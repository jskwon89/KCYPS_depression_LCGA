from __future__ import annotations

import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.cluster import KMeans
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.metrics import silhouette_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
FIG = OUT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)


ID_COLS = ["ID", "HID", "PID"]
WAVES = range(1, 8)


def strip_wave_suffix(col: str, wave: int) -> str:
    suffix = f"w{wave}"
    if col in ID_COLS:
        return col
    if col.endswith(suffix):
        return col[: -len(suffix)]
    return col


def read_wave(kind: str, wave: int) -> pd.DataFrame:
    path = ROOT / f"KCYPS2018m1{kind}w{wave}.dta"
    df = pd.read_stata(path, convert_categoricals=False)
    df = df.rename(columns={c: strip_wave_suffix(c, wave) for c in df.columns})
    df["wave"] = wave
    return df


def to_num(x: pd.DataFrame | pd.Series) -> pd.DataFrame | pd.Series:
    return x.apply(pd.to_numeric, errors="coerce") if isinstance(x, pd.DataFrame) else pd.to_numeric(x, errors="coerce")


def valid_items(df: pd.DataFrame, items: list[str], low: float, high: float) -> pd.DataFrame:
    existing = [c for c in items if c in df.columns]
    x = to_num(df[existing]).copy()
    return x.where(x.ge(low) & x.le(high))


def scale_mean(
    df: pd.DataFrame,
    items: list[str],
    low: float,
    high: float,
    reverse: set[str] | None = None,
    min_prop: float = 0.5,
) -> pd.Series:
    x = valid_items(df, items, low, high)
    reverse = reverse or set()
    for col in [c for c in x.columns if c in reverse]:
        x[col] = low + high - x[col]
    min_n = max(1, math.ceil(len(items) * min_prop))
    return x.mean(axis=1).where(x.notna().sum(axis=1).ge(min_n))


def alpha_from_items(
    df: pd.DataFrame,
    items: list[str],
    low: float,
    high: float,
    reverse: set[str] | None = None,
    min_prop: float = 0.8,
) -> float:
    x = valid_items(df, items, low, high)
    reverse = reverse or set()
    for col in [c for c in x.columns if c in reverse]:
        x[col] = low + high - x[col]
    min_n = max(2, math.ceil(len(items) * min_prop))
    x = x[x.notna().sum(axis=1).ge(min_n)]
    if x.shape[0] < 30 or x.shape[1] < 2:
        return np.nan
    x = x.fillna(x.mean())
    total = x.sum(axis=1)
    denom = total.var(ddof=1)
    if denom <= 0:
        return np.nan
    k = x.shape[1]
    return float(k / (k - 1) * (1 - x.var(ddof=1).sum() / denom))


def delinquency_scores(df: pd.DataFrame, items: list[str], min_prop: float = 0.5) -> tuple[pd.Series, pd.Series, pd.Series]:
    x = valid_items(df, items, 1, 6)
    min_n = max(1, math.ceil(len(items) * min_prop))
    usable = x.notna().sum(axis=1).ge(min_n)
    freq = (x - 1).mean(axis=1).where(usable)
    count = x.gt(1).sum(axis=1).where(usable)
    any_ = x.gt(1).any(axis=1).astype(float).where(usable)
    return freq, count, any_


def zscore(s: pd.Series) -> pd.Series:
    s = to_num(s)
    sd = s.std(ddof=0)
    if pd.isna(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean()) / sd


def weighted_mean(x: pd.Series, w: pd.Series | None = None) -> float:
    x = to_num(x)
    if w is None:
        return float(x.mean())
    w = to_num(w)
    m = x.notna() & w.notna() & (w > 0)
    if m.sum() == 0:
        return np.nan
    return float(np.average(x[m], weights=w[m]))


Y_ITEMS = {
    "offline": [f"YDLQ1A{i:02d}" for i in range(1, 16)],
    "cyber": [f"YDLQ2A{i:02d}" for i in range(1, 16)],
    "smartphone_dependence": [f"YMDA1C{i:02d}" for i in range(1, 16)],
    "parent_warmth": [f"YFAM2A{i:02d}" for i in range(1, 5)],
    "parent_rejection": [f"YFAM2B{i:02d}" for i in range(1, 5)],
    "parent_autonomy": [f"YFAM2C{i:02d}" for i in range(1, 5)],
    "parent_coercion": [f"YFAM2D{i:02d}" for i in range(1, 5)],
    "parent_structure": [f"YFAM2E{i:02d}" for i in range(1, 5)],
    "parent_inconsistency": [f"YFAM2F{i:02d}" for i in range(1, 5)],
    "positive_peer": [f"YEDU2A{i:02d}" for i in range(1, 9)],
    "negative_peer": [f"YEDU2A{i:02d}" for i in range(9, 14)],
    "teacher_relation": [f"YEDU3A{i:02d}" for i in range(1, 15)],
    "attention_problem": [f"YPSY4A{i:02d}" for i in range(1, 8)],
    "aggression": [f"YPSY4B{i:02d}" for i in range(1, 7)],
    "depression": [f"YPSY4E{i:02d}" for i in range(1, 11)],
    "self_esteem": [f"YPSY3A{i:02d}" for i in range(1, 11)],
    "grit": [f"YPSY7A{i:02d}" for i in range(1, 9)],
    "academic_engagement": [f"YINT2A{i:02d}" for i in range(1, 17)],
    "academic_helplessness": [f"YINT2B{i:02d}" for i in range(1, 17)],
    "life_satisfaction": [f"YPSY1A{i:02d}" for i in range(1, 6)],
}

REVERSE = {
    "smartphone_dependence": {"YMDA1C05", "YMDA1C10", "YMDA1C15"},
    "self_esteem": {"YPSY3A02", "YPSY3A05", "YPSY3A06", "YPSY3A08", "YPSY3A09"},
    "grit": {"YPSY7A01", "YPSY7A03", "YPSY7A05", "YPSY7A06"},
}

CORE_CONSTRUCTS = [
    "cyber_freq",
    "offline_freq",
    "smartphone_dependence",
    "parent_warmth",
    "parent_rejection",
    "parent_coercion",
    "parent_inconsistency",
    "positive_peer",
    "negative_peer",
    "teacher_relation",
    "attention_problem",
    "aggression",
    "depression",
    "self_esteem",
    "grit",
    "academic_engagement",
    "academic_helplessness",
    "life_satisfaction",
]


def build_panel() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    reliabilities = []

    for wave in WAVES:
        y = read_wave("Y", wave)
        p = read_wave("P", wave)
        p_keep = [
            c
            for c in [
                "ID",
                "PINCOME",
                "PINCOME1",
                "PSCHOOL1",
                "PSCHOOL2",
                "PWORKW1",
                "PWORKW2",
            ]
            if c in p.columns
        ]
        df = y.merge(p[p_keep], on="ID", how="left", suffixes=("", "_p"))

        out = pd.DataFrame({"ID": df["ID"], "HID": df["HID"], "PID": df["PID"], "wave": wave})
        out["survey_youth"] = to_num(df["SURVEY1"]) if "SURVEY1" in df.columns else 1
        out["survey_parent"] = to_num(df["SURVEY2"]) if "SURVEY2" in df.columns else np.nan
        out["gender_female"] = to_num(df["YGENDER"]).map({1.0: 0.0, 2.0: 1.0})
        out["urban_size"] = to_num(df["ARA2A"]).where(to_num(df["ARA2A"]).between(1, 3))
        out["school_type"] = to_num(df["SCLTYP"]).where(to_num(df["SCLTYP"]).between(1, 7)) if "SCLTYP" in df.columns else np.nan
        out["weight_cross_std"] = to_num(df["WEIGHTA2"]) if "WEIGHTA2" in df.columns else np.nan
        out["weight_long_std"] = to_num(df["WEIGHTB2"]) if "WEIGHTB2" in df.columns else np.nan

        cyber_freq, cyber_count, cyber_any = delinquency_scores(df, Y_ITEMS["cyber"])
        offline_freq, offline_count, offline_any = delinquency_scores(df, Y_ITEMS["offline"])
        out["cyber_freq"] = cyber_freq
        out["cyber_count"] = cyber_count
        out["cyber_any"] = cyber_any
        out["offline_freq"] = offline_freq
        out["offline_count"] = offline_count
        out["offline_any"] = offline_any

        for name, items in Y_ITEMS.items():
            if name in {"cyber", "offline"}:
                reliabilities.append(
                    {
                        "wave": wave,
                        "scale": f"{name}_frequency",
                        "alpha": alpha_from_items(df, items, 1, 6),
                    }
                )
                continue
            out[name] = scale_mean(df, items, 1, 4, reverse=REVERSE.get(name, set()))
            reliabilities.append(
                {
                    "wave": wave,
                    "scale": name,
                    "alpha": alpha_from_items(df, items, 1, 4, reverse=REVERSE.get(name, set())),
                }
            )

        out["academic_achievement"] = to_num(df["YINT1A00"]).where(to_num(df["YINT1A00"]).between(1, 5))
        out["school_satisfaction"] = to_num(df["YEDU1A00"]).where(to_num(df["YEDU1A00"]).between(1, 5)) if "YEDU1A00" in df.columns else np.nan
        out["parent_income"] = to_num(df["PINCOME"]).where(to_num(df["PINCOME"]).between(1, 12))
        out["subjective_ses"] = to_num(df["PINCOME1"]).where(to_num(df["PINCOME1"]).between(1, 5))
        father_ed = to_num(df["PSCHOOL1"]).where(to_num(df["PSCHOOL1"]).between(1, 7))
        mother_ed = to_num(df["PSCHOOL2"]).where(to_num(df["PSCHOOL2"]).between(1, 7))
        out["parent_education_max"] = pd.concat([father_ed, mother_ed], axis=1).max(axis=1)
        out["father_work"] = to_num(df["PWORKW1"]).map({1.0: 1.0, 2.0: 0.0})
        out["mother_work"] = to_num(df["PWORKW2"]).map({1.0: 1.0, 2.0: 0.0})
        rows.append(out)

    panel = pd.concat(rows, ignore_index=True).sort_values(["ID", "wave"])
    reliabilities_df = pd.DataFrame(reliabilities)
    return panel, reliabilities_df


def make_descriptives(panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for wave, g in panel.groupby("wave"):
        rows.append(
            {
                "wave": wave,
                "n_rows": len(g),
                "n_cyber_valid": int(g["cyber_freq"].notna().sum()),
                "n_offline_valid": int(g["offline_freq"].notna().sum()),
                "cyber_any": g["cyber_any"].mean(),
                "offline_any": g["offline_any"].mean(),
                "cyber_freq_mean": g["cyber_freq"].mean(),
                "offline_freq_mean": g["offline_freq"].mean(),
                "smartphone_dependence": g["smartphone_dependence"].mean(),
                "depression": g["depression"].mean(),
                "aggression": g["aggression"].mean(),
                "parent_rejection": g["parent_rejection"].mean(),
                "negative_peer": g["negative_peer"].mean(),
                "weighted_cyber_any": weighted_mean(g["cyber_any"], g["weight_cross_std"]),
                "weighted_offline_any": weighted_mean(g["offline_any"], g["weight_cross_std"]),
            }
        )
    return pd.DataFrame(rows)


def make_transition_data(panel: pd.DataFrame) -> pd.DataFrame:
    panel = panel.sort_values(["ID", "wave"]).copy()
    lag_vars = CORE_CONSTRUCTS + [
        "cyber_any",
        "offline_any",
        "cyber_count",
        "offline_count",
        "academic_achievement",
        "school_satisfaction",
        "parent_income",
        "subjective_ses",
        "parent_education_max",
        "gender_female",
        "urban_size",
        "school_type",
    ]
    for v in lag_vars:
        if v in panel.columns:
            panel[f"lag_{v}"] = panel.groupby("ID")[v].shift(1)
    return panel[panel["wave"].ge(2)].copy()


def within_fe_model(df: pd.DataFrame, outcome: str, predictors: list[str]) -> pd.DataFrame:
    cols = ["ID", "wave", outcome] + predictors
    d = df[cols].dropna().copy()
    counts = d.groupby("ID")[outcome].transform("count")
    d = d[counts.ge(2)].copy()

    y = zscore(d[outcome])
    x = pd.DataFrame({p: zscore(d[p]) for p in predictors}, index=d.index)
    wave_dummies = pd.get_dummies(d["wave"].astype(int), prefix="wave", drop_first=True, dtype=float)
    x = pd.concat([x, wave_dummies], axis=1)

    group = d["ID"]
    y_w = y - y.groupby(group).transform("mean")
    x_w = x - x.groupby(group).transform("mean")

    keep = y_w.notna() & x_w.notna().all(axis=1)
    y_w = y_w.loc[keep]
    x_w = x_w.loc[keep]
    group = group.loc[keep]
    model = sm.OLS(y_w, x_w).fit(cov_type="cluster", cov_kwds={"groups": group})

    table = []
    for name in predictors:
        table.append(
            {
                "outcome": outcome,
                "term": name,
                "coef_std": model.params.get(name, np.nan),
                "se_cluster": model.bse.get(name, np.nan),
                "p": model.pvalues.get(name, np.nan),
                "n_obs": int(model.nobs),
                "n_ids": int(group.nunique()),
                "r2_within": float(model.rsquared),
            }
        )
    return pd.DataFrame(table)


def trajectory_clustering(panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    wide = panel.pivot(index="ID", columns="wave", values="cyber_freq").sort_index(axis=1)
    keep = wide.notna().sum(axis=1).ge(5)
    wide = wide[keep]
    imp = wide.T.interpolate(limit_direction="both").T
    imp = imp.apply(lambda row: row.fillna(row.mean()), axis=1)
    imp = imp.fillna(0)
    x = StandardScaler().fit_transform(imp)

    fit_rows = []
    labels_by_k = {}
    max_k = 6
    for k in range(2, max_k + 1):
        km = KMeans(n_clusters=k, random_state=20260507, n_init=50)
        labels = km.fit_predict(x)
        sil = silhouette_score(x, labels) if len(set(labels)) > 1 else np.nan
        fit_rows.append({"k": k, "silhouette": sil, "inertia": km.inertia_})
        labels_by_k[k] = labels

    fit = pd.DataFrame(fit_rows)
    best_k = int(fit.sort_values(["silhouette", "k"], ascending=[False, True]).iloc[0]["k"])
    labels = labels_by_k[best_k]
    assign = pd.DataFrame({"ID": imp.index, "cluster": labels})
    traj = imp.copy()
    traj["cluster"] = labels
    means = traj.groupby("cluster").mean().reset_index()

    label_map = {}
    for _, row in means.iterrows():
        vals = row[[1, 2, 3, 4, 5, 6, 7]].astype(float).to_numpy()
        avg = vals.mean()
        slope = vals[-1] - vals[0]
        peak = vals.max()
        if avg < 0.03 and peak < 0.08:
            name = "low-stable"
        elif slope > 0.08:
            name = "escalating"
        elif slope < -0.08:
            name = "declining"
        elif avg >= 0.12:
            name = "persistent-high"
        else:
            name = "intermittent/moderate"
        label_map[int(row["cluster"])] = name
    assign["cluster_label"] = assign["cluster"].map(label_map)
    means.insert(1, "cluster_label", means["cluster"].map(label_map))
    sizes = assign["cluster_label"].value_counts().rename_axis("cluster_label").reset_index(name="n")
    sizes["pct"] = sizes["n"] / sizes["n"].sum()
    means = means.merge(sizes, on="cluster_label", how="left")
    return fit, assign, means


def ml_prediction(transitions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    target = "cyber_any"
    feature_vars = [
        "lag_cyber_any",
        "lag_cyber_freq",
        "lag_offline_any",
        "lag_offline_freq",
        "lag_smartphone_dependence",
        "lag_depression",
        "lag_aggression",
        "lag_attention_problem",
        "lag_parent_warmth",
        "lag_parent_rejection",
        "lag_parent_coercion",
        "lag_parent_inconsistency",
        "lag_positive_peer",
        "lag_negative_peer",
        "lag_teacher_relation",
        "lag_self_esteem",
        "lag_grit",
        "lag_academic_engagement",
        "lag_academic_helplessness",
        "lag_life_satisfaction",
        "lag_academic_achievement",
        "lag_school_satisfaction",
        "lag_parent_income",
        "lag_subjective_ses",
        "lag_parent_education_max",
        "lag_gender_female",
        "lag_urban_size",
        "lag_school_type",
    ]
    d = transitions[["ID", "wave", target] + feature_vars].copy()
    d = d[d[target].notna()].copy()
    train = d[d["wave"].between(2, 6)].copy()
    test = d[d["wave"].eq(7)].copy()
    x_train, y_train = train[feature_vars], train[target].astype(int)
    x_test, y_test = test[feature_vars], test[target].astype(int)

    models = {
        "lag-only": Pipeline(
            [
                ("impute", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(max_iter=2000, class_weight="balanced", random_state=20260507),
                ),
            ]
        ),
        "lasso-logit": Pipeline(
            [
                ("impute", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=4000,
                        class_weight="balanced",
                        penalty="l1",
                        solver="saga",
                        C=0.15,
                        random_state=20260507,
                    ),
                ),
            ]
        ),
        "random-forest": Pipeline(
            [
                ("impute", SimpleImputer(strategy="median")),
                (
                    "model",
                    RandomForestClassifier(
                        n_estimators=700,
                        min_samples_leaf=8,
                        class_weight="balanced_subsample",
                        random_state=20260507,
                        n_jobs=-1,
                    ),
                ),
            ]
        ),
        "hist-gradient-boosting": Pipeline(
            [
                ("impute", SimpleImputer(strategy="median")),
                (
                    "model",
                    HistGradientBoostingClassifier(
                        max_iter=300,
                        learning_rate=0.04,
                        l2_regularization=0.04,
                        max_leaf_nodes=15,
                        random_state=20260507,
                    ),
                ),
            ]
        ),
    }

    # The lag-only benchmark uses only prior cyber/offline variables.
    lag_features = ["lag_cyber_any", "lag_cyber_freq", "lag_offline_any", "lag_offline_freq"]

    rows = []
    fitted = {}
    for name, model in models.items():
        feats = lag_features if name == "lag-only" else feature_vars
        model.fit(x_train[feats], y_train)
        prob = model.predict_proba(x_test[feats])[:, 1]
        rows.append(
            {
                "model": name,
                "n_train": len(train),
                "n_test": len(test),
                "test_event_rate": y_test.mean(),
                "roc_auc": roc_auc_score(y_test, prob),
                "average_precision": average_precision_score(y_test, prob),
                "brier": brier_score_loss(y_test, prob),
            }
        )
        fitted[name] = (model, feats)

    metrics = pd.DataFrame(rows).sort_values("roc_auc", ascending=False)
    best_name = metrics.iloc[0]["model"]
    best_model, best_feats = fitted[best_name]
    perm = permutation_importance(
        best_model,
        x_test[best_feats],
        y_test,
        scoring="roc_auc",
        n_repeats=20,
        random_state=20260507,
        n_jobs=-1,
    )
    importance = (
        pd.DataFrame(
            {
                "model": best_name,
                "feature": best_feats,
                "importance_mean": perm.importances_mean,
                "importance_sd": perm.importances_std,
            }
        )
        .sort_values("importance_mean", ascending=False)
        .reset_index(drop=True)
    )
    return metrics, importance


def gee_cyber_onset_model(transitions: pd.DataFrame) -> pd.DataFrame:
    """Discrete-time risk model among youth without cyber delinquency in t-1."""
    risk = transitions[(transitions["lag_cyber_any"].eq(0)) & transitions["cyber_any"].notna()].copy()
    predictors_cont = [
        "lag_smartphone_dependence",
        "lag_depression",
        "lag_aggression",
        "lag_attention_problem",
        "lag_parent_warmth",
        "lag_parent_rejection",
        "lag_parent_coercion",
        "lag_parent_inconsistency",
        "lag_negative_peer",
        "lag_teacher_relation",
        "lag_self_esteem",
        "lag_grit",
        "lag_academic_engagement",
        "lag_academic_helplessness",
        "lag_life_satisfaction",
        "lag_offline_freq",
        "lag_parent_income",
        "lag_parent_education_max",
    ]
    predictors_bin = ["lag_gender_female"]
    cols = ["ID", "wave", "cyber_any"] + predictors_cont + predictors_bin
    d = risk[cols].dropna().copy()
    for c in predictors_cont:
        d[c] = zscore(d[c])
    x = pd.concat(
        [
            d[predictors_cont + predictors_bin],
            pd.get_dummies(d["wave"].astype(int), prefix="wave", drop_first=True, dtype=float),
        ],
        axis=1,
    )
    x = sm.add_constant(x)
    y = d["cyber_any"].astype(int)
    model = sm.GEE(
        y,
        x,
        groups=d["ID"],
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Exchangeable(),
    )
    res = model.fit()
    rows = []
    for term in predictors_cont + predictors_bin:
        b = res.params[term]
        se = res.bse[term]
        rows.append(
            {
                "term": term,
                "coef_logit": b,
                "se": se,
                "OR": np.exp(b),
                "ci_low": np.exp(b - 1.96 * se),
                "ci_high": np.exp(b + 1.96 * se),
                "p": res.pvalues[term],
                "n_obs": int(len(d)),
                "n_ids": int(d["ID"].nunique()),
                "event_rate": float(y.mean()),
                "working_corr": float(np.atleast_1d(res.cov_struct.dep_params)[0]),
            }
        )
    return pd.DataFrame(rows).sort_values("p")


def baseline_cluster_prediction(panel: pd.DataFrame, assignments: pd.DataFrame) -> pd.DataFrame:
    base = panel[panel["wave"].eq(1)].merge(assignments[["ID", "cluster_label"]], on="ID", how="inner")
    predictors = [
        "cyber_any",
        "offline_any",
        "smartphone_dependence",
        "depression",
        "aggression",
        "attention_problem",
        "parent_warmth",
        "parent_rejection",
        "parent_inconsistency",
        "negative_peer",
        "teacher_relation",
        "self_esteem",
        "grit",
        "academic_engagement",
        "parent_income",
        "parent_education_max",
        "gender_female",
    ]
    rows = []
    for group, g in base.groupby("cluster_label"):
        rec = {"cluster_label": group, "n": len(g)}
        for p in predictors:
            rec[p] = g[p].mean()
        rows.append(rec)
    return pd.DataFrame(rows).sort_values("n", ascending=False)


def save_figures(desc: pd.DataFrame, traj_means: pd.DataFrame, ml_imp: pd.DataFrame) -> None:
    plt.figure(figsize=(8, 4.8))
    plt.plot(desc["wave"], desc["cyber_any"], marker="o", label="Any cyber delinquency")
    plt.plot(desc["wave"], desc["offline_any"], marker="o", label="Any offline delinquency")
    plt.xlabel("Wave")
    plt.ylabel("Prevalence")
    plt.ylim(0, max(desc["offline_any"].max(), desc["cyber_any"].max()) * 1.25)
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIG / "trend_prevalence.png", dpi=220)
    plt.close()

    plt.figure(figsize=(8, 4.8))
    for _, row in traj_means.iterrows():
        vals = [row[w] for w in WAVES]
        label = f"{row['cluster_label']} (n={int(row['n'])})"
        plt.plot(list(WAVES), vals, marker="o", label=label)
    plt.xlabel("Wave")
    plt.ylabel("Mean cyber delinquency frequency (0-5)")
    plt.grid(alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG / "cyber_trajectory_clusters.png", dpi=220)
    plt.close()

    top = ml_imp.head(12).iloc[::-1]
    plt.figure(figsize=(8, 5.2))
    plt.barh(top["feature"], top["importance_mean"], xerr=top["importance_sd"], color="#4C78A8")
    plt.xlabel("Permutation importance in ROC AUC")
    plt.tight_layout()
    plt.savefig(FIG / "ml_permutation_importance.png", dpi=220)
    plt.close()


def fmt_p(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "<.001"
    return f"{p:.3f}"


def write_report(
    desc: pd.DataFrame,
    reliability: pd.DataFrame,
    fe: pd.DataFrame,
    fe_no_prior: pd.DataFrame,
    gee_onset: pd.DataFrame,
    traj_fit: pd.DataFrame,
    traj_means: pd.DataFrame,
    ml_metrics: pd.DataFrame,
    ml_imp: pd.DataFrame,
    cluster_base: pd.DataFrame,
) -> None:
    cyber_fe = fe[fe["outcome"].eq("cyber_freq")]
    cyber_top = cyber_fe.reindex(cyber_fe["coef_std"].abs().sort_values(ascending=False).index).head(8)
    offline_fe = fe[fe["outcome"].eq("offline_freq")]
    offline_top = offline_fe.reindex(offline_fe["coef_std"].abs().sort_values(ascending=False).index).head(8)
    cyber_no_prior = fe_no_prior[fe_no_prior["outcome"].eq("cyber_freq")]
    cyber_no_prior_top = cyber_no_prior.sort_values("p").head(8)
    onset_top = gee_onset.head(10)

    alpha_summary = (
        reliability.groupby("scale")["alpha"].agg(["mean", "min", "max"]).reset_index().sort_values("mean", ascending=False)
    )

    lines: list[str] = []
    lines.append("# KCYPS 2018 중1 코호트 1~7차 패널 분석: SSCI급 논문 설계와 예비결과")
    lines.append("")
    lines.append("## 1. 데이터와 설계")
    lines.append(
        "- 자료: `KCYPS 2018` 중1 코호트 원패널 청소년(Y) 및 보호자(P) 1~7차 DTA, 유저가이드 2025.12판, 코드북 Excel."
    )
    lines.append(
        "- 분석 단위: 청소년 개인-파형 패널. 청소년 파일은 각 파형 2,590행으로 구성되어 있으며, 문항 비응답/미참여는 구성척도 산출에서 제외했다."
    )
    lines.append(
        "- 유저가이드 기준 표집틀은 2018년 초4·중1 재학생이며, 지역과 도시규모를 층화축으로 한 다단계층화집락추출 설계다. 중1 코호트 원패널은 2,590명이고 7차 조사 완료자는 2,073명, 유지율은 80.0%다."
    )
    lines.append(
        "- 코드에는 횡단면 표준화 가중치와 종단면 표준화 가중치가 포함되어 있어, 본문 결과는 비가중 예비결과로 제시하고 기술통계에는 가중 prevalence를 함께 산출했다."
    )
    lines.append(
        "- 핵심 결과변수: 현실비행 15문항과 사이버비행 15문항을 각각 0~5 빈도 지수와 경험 여부로 재코딩했다."
    )
    lines.append(
        "- 핵심 설명변수: 스마트폰 의존, 부모 양육태도, 또래/교사관계, 정서문제, 자아존중감, 그릿, 학업태도, SES."
    )
    lines.append("")

    lines.append("## 2. 가장 강한 논문 아이디어")
    lines.append(
        "**Adolescent cyber delinquency as a developmental cascade:** 스마트폰 의존과 정서·관계 요인이 현실비행과 독립적으로 사이버비행의 청소년기 변화와 이행을 예측하는가?"
    )
    lines.append(
        "이 주제가 유망한 이유는 7개 파형 전체에 사이버비행 15문항이 반복 측정되어 있고, 스마트폰 의존·정서문제·부모/또래/교사 관계도 같은 기간 반복 측정되어 패널 고정효과, 교차지연 예측, 궤적군집, 머신러닝 예측을 한 논문 안에 결합할 수 있기 때문이다."
    )
    lines.append("")

    lines.append("## 3. 기술통계 핵심")
    lines.append(desc.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")
    lines.append("신뢰도(Cronbach alpha) 요약:")
    lines.append(alpha_summary.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")

    lines.append("## 4. 패널 고정효과 결과")
    lines.append(
        "개인 고정효과와 파형 고정효과를 두고, t-1 시점의 표준화된 예측변수가 t 시점의 표준화된 비행 빈도를 예측하는 모형을 추정했다. 계수는 개인 내 변화 기준의 표준화 효과다. 짧은 패널에서 lagged outcome을 넣은 동적 FE는 평균회귀가 강하게 보일 수 있어, 직전 비행을 제외한 모형과 함께 해석한다."
    )
    lines.append("")
    lines.append("직전 비행을 제외한 사이버비행 FE 모형:")
    lines.append(
        cyber_no_prior_top.assign(p=cyber_no_prior_top["p"].map(fmt_p))[
            ["term", "coef_std", "se_cluster", "p", "n_obs", "n_ids", "r2_within"]
        ].to_markdown(index=False, floatfmt=".3f")
    )
    lines.append("")
    lines.append("사이버비행 상위 예측요인:")
    lines.append(
        cyber_top.assign(p=cyber_top["p"].map(fmt_p))[
            ["term", "coef_std", "se_cluster", "p", "n_obs", "n_ids", "r2_within"]
        ].to_markdown(index=False, floatfmt=".3f")
    )
    lines.append("")
    lines.append("현실비행 상위 예측요인:")
    lines.append(
        offline_top.assign(p=offline_top["p"].map(fmt_p))[
            ["term", "coef_std", "se_cluster", "p", "n_obs", "n_ids", "r2_within"]
        ].to_markdown(index=False, floatfmt=".3f")
    )
    lines.append("")

    lines.append("## 5. 사이버비행 진입 위험모형")
    lines.append(
        "직전 파형에 사이버비행 경험이 없던 청소년만 위험집합으로 두고, 다음 파형의 새 사이버비행 경험 여부를 GEE logit으로 추정했다. 연속 예측변수는 표준화했으므로 OR은 1 SD 증가의 효과이며, 성별은 여학생=1이다."
    )
    onset_display = onset_top.copy()
    onset_display["p"] = onset_display["p"].map(fmt_p)
    lines.append(
        onset_display[["term", "OR", "ci_low", "ci_high", "p", "n_obs", "n_ids", "event_rate"]].to_markdown(
            index=False, floatfmt=".3f"
        )
    )
    lines.append(
        "해석: 남학생(여학생 OR<1), 높은 스마트폰 의존, 공격성, 직전 현실비행, 부모 거부가 사이버비행 진입 위험을 높이는 가장 일관적인 신호다."
    )
    lines.append("")

    lines.append("## 6. 궤적군집")
    lines.append("K-means 기반 반복측정 궤적군집을 2~6개 군으로 비교했고 silhouette를 기준으로 예비 군집 수를 선택했다.")
    lines.append(traj_fit.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")
    lines.append("선택 군집별 사이버비행 평균 궤적:")
    lines.append(traj_means.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")
    lines.append("군집별 1차년도 위험요인 평균:")
    lines.append(cluster_base.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")

    lines.append("## 7. 머신러닝 예측")
    lines.append(
        "1~5차 전이자료로 학습하고 6→7차 전이를 holdout test로 두어, 7차 사이버비행 경험 여부를 예측했다. 비교 지표는 ROC-AUC, average precision, Brier score다."
    )
    lines.append(ml_metrics.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")
    lines.append("최상위 모형의 permutation importance:")
    lines.append(ml_imp.head(15).to_markdown(index=False, floatfmt=".4f"))
    lines.append("")

    lines.append("## 8. 투고용 분석 설계")
    lines.append("1. 측정: 15문항 사이버비행 지수와 현실비행 지수의 파형별 신뢰도와 불변성 점검.")
    lines.append("2. 기술 분석: 성별/SES/학교유형별 추세, 가중치 적용 민감도.")
    lines.append("3. 원인 추론에 가까운 패널 분석: 개인 고정효과 + 파형 고정효과 + lagged predictors.")
    lines.append("4. 발달 이질성: LCGA/GMM 또는 본 예비분석의 군집을 바탕으로 잠재궤적모형 확정.")
    lines.append("5. 예측 타당도: 시간차 holdout ML로 결과의 외적 예측 성능 제시.")
    lines.append("6. 강건성: 현실비행 통제, 성별 상호작용, zero-inflated/negative-binomial 대체모형, 완전사례와 다중대체 비교.")
    lines.append("")
    lines.append("## 9. 그림")
    lines.append("![trend](figures/trend_prevalence.png)")
    lines.append("![trajectory](figures/cyber_trajectory_clusters.png)")
    lines.append("![importance](figures/ml_permutation_importance.png)")
    lines.append("")

    (OUT / "kcyps_ssci_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    panel, reliability = build_panel()
    panel.to_csv(OUT / "analytic_panel.csv", index=False, encoding="utf-8-sig")
    reliability.to_csv(OUT / "scale_reliability.csv", index=False, encoding="utf-8-sig")

    desc = make_descriptives(panel)
    desc.to_csv(OUT / "descriptives_by_wave.csv", index=False, encoding="utf-8-sig")

    transitions = make_transition_data(panel)
    transitions.to_csv(OUT / "transition_panel_lagged.csv", index=False, encoding="utf-8-sig")

    predictors = [
        "lag_offline_freq",
        "lag_cyber_freq",
        "lag_smartphone_dependence",
        "lag_depression",
        "lag_aggression",
        "lag_attention_problem",
        "lag_parent_warmth",
        "lag_parent_rejection",
        "lag_parent_coercion",
        "lag_parent_inconsistency",
        "lag_negative_peer",
        "lag_teacher_relation",
        "lag_self_esteem",
        "lag_grit",
        "lag_academic_engagement",
        "lag_academic_helplessness",
        "lag_life_satisfaction",
        "lag_parent_income",
        "lag_parent_education_max",
    ]
    fe_tables = [
        within_fe_model(transitions, "cyber_freq", predictors),
        within_fe_model(transitions, "offline_freq", predictors),
    ]
    fe = pd.concat(fe_tables, ignore_index=True)
    fe.to_csv(OUT / "fixed_effects_lagged_models.csv", index=False, encoding="utf-8-sig")

    no_prior_predictors = [
        p for p in predictors if p not in {"lag_offline_freq", "lag_cyber_freq"}
    ]
    fe_no_prior = pd.concat(
        [
            within_fe_model(transitions, "cyber_freq", no_prior_predictors),
            within_fe_model(transitions, "offline_freq", no_prior_predictors),
        ],
        ignore_index=True,
    )
    fe_no_prior.to_csv(OUT / "fixed_effects_no_prior_delinquency.csv", index=False, encoding="utf-8-sig")

    gee_onset = gee_cyber_onset_model(transitions)
    gee_onset.to_csv(OUT / "gee_cyber_onset_model.csv", index=False, encoding="utf-8-sig")

    traj_fit, assignments, traj_means = trajectory_clustering(panel)
    traj_fit.to_csv(OUT / "trajectory_fit.csv", index=False, encoding="utf-8-sig")
    assignments.to_csv(OUT / "trajectory_assignments.csv", index=False, encoding="utf-8-sig")
    traj_means.to_csv(OUT / "trajectory_means.csv", index=False, encoding="utf-8-sig")
    cluster_base = baseline_cluster_prediction(panel, assignments)
    cluster_base.to_csv(OUT / "trajectory_baseline_profiles.csv", index=False, encoding="utf-8-sig")

    ml_metrics, ml_imp = ml_prediction(transitions)
    ml_metrics.to_csv(OUT / "ml_prediction_metrics.csv", index=False, encoding="utf-8-sig")
    ml_imp.to_csv(OUT / "ml_permutation_importance.csv", index=False, encoding="utf-8-sig")

    save_figures(desc, traj_means, ml_imp)
    write_report(desc, reliability, fe, fe_no_prior, gee_onset, traj_fit, traj_means, ml_metrics, ml_imp, cluster_base)

    summary = {
        "n_panel_rows": int(len(panel)),
        "n_ids": int(panel["ID"].nunique()),
        "outputs": [
            str(OUT / "analytic_panel.csv"),
            str(OUT / "kcyps_ssci_report.md"),
            str(OUT / "fixed_effects_lagged_models.csv"),
            str(OUT / "ml_prediction_metrics.csv"),
            str(OUT / "trajectory_means.csv"),
        ],
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
