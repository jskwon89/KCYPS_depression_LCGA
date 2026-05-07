from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CODEBOOK = ROOT / "KCYPS 2018_Codebook.xlsx"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)


def load_layout(sheet: str = "LAYOUT_m") -> pd.DataFrame:
    df = pd.read_excel(CODEBOOK, sheet_name=sheet, header=None, skiprows=5)
    df = df.iloc[:, 1:15].copy()
    df.columns = [
        "domain",
        "mid",
        "sub",
        "item",
        "factor",
        "desc",
        "var",
        "w1",
        "w2",
        "w3",
        "w4",
        "w5",
        "w6",
        "w7",
    ]
    df = df[df["var"].notna()].copy()
    for col in df.columns:
        df[col] = df[col].where(df[col].notna(), "")
    df["waves_present"] = df[[f"w{i}" for i in range(1, 8)]].apply(
        lambda row: [i for i, v in enumerate(row, start=1) if str(v).strip() != ""],
        axis=1,
    )
    return df


def compact_records(df: pd.DataFrame, n: int = 200) -> list[dict[str, object]]:
    cols = ["domain", "mid", "sub", "item", "factor", "desc", "var", "waves_present"]
    return df[cols].head(n).to_dict(orient="records")


def main() -> None:
    layout = load_layout()
    layout.to_csv(OUT / "layout_m_parsed.csv", index=False, encoding="utf-8-sig")

    summary = {
        "n_layout_rows": int(len(layout)),
        "domain_counts": layout["domain"].value_counts().to_dict(),
        "mid_counts": layout["mid"].value_counts().to_dict(),
        "sub_counts_top50": layout["sub"].value_counts().head(50).to_dict(),
    }

    focus_terms = [
        "비행",
        "공격",
        "우울",
        "불안",
        "스마트",
        "휴대",
        "인터넷",
        "자아",
        "부모",
        "방임",
        "학대",
        "친구",
        "교사",
        "학교",
        "진로",
        "학업",
        "자살",
        "삶",
        "행복",
        "공동체",
        "그릿",
        "성취",
        "수면",
        "운동",
        "주의",
        "피해",
        "가해",
    ]
    focus = {}
    searchable_cols = ["domain", "mid", "sub", "item", "factor", "desc", "var"]
    haystack = layout[searchable_cols].astype(str).agg(" ".join, axis=1)
    for term in focus_terms:
        rows = layout[haystack.str.contains(term, regex=False, na=False)]
        focus[term] = {
            "n": int(len(rows)),
            "records": compact_records(rows, n=80),
        }

    by_mid = {}
    for mid in layout["mid"].drop_duplicates():
        if mid and mid != "-":
            rows = layout[layout["mid"].eq(mid)]
            by_mid[mid] = compact_records(rows, n=120)

    result = {"summary": summary, "focus": focus, "by_mid": by_mid}
    with open(OUT / "kcyps_codebook_exploration.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    for term in focus_terms:
        print(f"\n## {term} ({focus[term]['n']})")
        for rec in focus[term]["records"][:20]:
            print(
                f"{rec['var']}: {rec['mid']} / {rec['sub']} / "
                f"{rec['item']} / {rec['factor']} / {rec['desc']} / {rec['waves_present']}"
            )


if __name__ == "__main__":
    main()
