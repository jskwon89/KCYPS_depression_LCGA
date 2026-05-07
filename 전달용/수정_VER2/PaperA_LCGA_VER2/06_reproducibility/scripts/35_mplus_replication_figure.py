# -*- coding: utf-8 -*-
"""
Create a supplementary figure summarizing Python vs Mplus LCGA replication.

The project root is inferred from this script's location so the script remains
portable when the delivery package is copied.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def locate_mplus_paths() -> tuple[Path, Path]:
    """Return input/output paths in either the analysis tree or delivery tree."""
    here = Path(__file__).resolve()
    for candidate in [here.parent, *here.parents]:
        analysis_mplus = candidate / "output" / "paperA_lcga" / "primary_9item" / "mplus"
        if (analysis_mplus / "mplus_vs_python_fit_summary_9item.csv").exists():
            return analysis_mplus, analysis_mplus.parent / "figures_mplus"

        delivery_mplus = candidate / "04_results_data" / "mplus"
        if (delivery_mplus / "mplus_vs_python_fit_summary_9item.csv").exists():
            return delivery_mplus, candidate / "03_figures" / "Mplus_validation"

    raise FileNotFoundError(
        "Could not locate Mplus replication CSVs in either output/.../mplus "
        "or 04_results_data/mplus."
    )


MPLUS, OUT = locate_mplus_paths()
OUT.mkdir(parents=True, exist_ok=True)


def main() -> None:
    fit = pd.read_csv(MPLUS / "mplus_vs_python_fit_summary_9item.csv")
    align = pd.read_csv(MPLUS / "python_vs_mplus_class_comparison_9item.csv")

    fit["K"] = fit["K"].astype(int)
    for col in ["BIC", "AIC", "aBIC"]:
        fit[f"{col}_absdiff"] = (fit[f"{col}_Mplus"] - fit[f"{col}_Python"]).abs()

    colors = {
        "Python EM": "#1f77b4",
        "Mplus 7": "#d62728",
        "BIC": "#1f77b4",
        "AIC": "#2ca02c",
        "aBIC": "#ff7f0e",
    }

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })

    fig = plt.figure(figsize=(14.5, 4.9), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.15])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # Panel A: BIC overlap.
    ax1.plot(
        fit["K"],
        fit["BIC_Python"],
        marker="o",
        lw=2.0,
        color=colors["Python EM"],
        label="Python EM",
    )
    ax1.plot(
        fit["K"],
        fit["BIC_Mplus"],
        marker="s",
        lw=1.6,
        color=colors["Mplus 7"],
        linestyle="--",
        label="Mplus 7",
    )
    ax1.axvline(4, color="0.45", linestyle=":", lw=1.1)
    ax1.annotate(
        "K=4 primary",
        xy=(4, fit.loc[fit["K"].eq(4), "BIC_Mplus"].iloc[0]),
        xytext=(4.25, fit["BIC_Mplus"].min() + 500),
        arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.8},
        fontsize=8,
        color="0.25",
    )
    ax1.set_title("A. BIC by number of classes", loc="left", fontweight="bold")
    ax1.set_xlabel("Number of latent classes (K)")
    ax1.set_ylabel("BIC")
    ax1.set_xticks(fit["K"])
    ax1.grid(alpha=0.25)
    ax1.legend(frameon=False, loc="upper right")

    # Panel B: absolute differences across all major information criteria.
    diff_long = fit.melt(
        id_vars="K",
        value_vars=["BIC_absdiff", "AIC_absdiff", "aBIC_absdiff"],
        var_name="criterion",
        value_name="abs_difference",
    )
    diff_long["criterion"] = diff_long["criterion"].str.replace("_absdiff", "", regex=False)
    offsets = {"AIC": -0.18, "BIC": 0.0, "aBIC": 0.18}
    for crit in ["AIC", "BIC", "aBIC"]:
        sub = diff_long[diff_long["criterion"] == crit]
        ax2.bar(
            sub["K"] + offsets[crit],
            sub["abs_difference"],
            width=0.16,
            color=colors[crit],
            label=crit,
        )
    ax2.set_title("B. |Mplus - Python| fit-index difference", loc="left", fontweight="bold")
    ax2.set_xlabel("K")
    ax2.set_ylabel("Absolute difference")
    ax2.set_xticks(fit["K"])
    ax2.set_ylim(0, max(0.003, diff_long["abs_difference"].max() * 1.25))
    ax2.grid(axis="y", alpha=0.25)
    ax2.legend(frameon=False, loc="upper right")

    # Panel C: class counts after label alignment.
    labels = [
        "C1\nStable-Low",
        "C2\nPersistent-High",
        "C3\nDecreasing",
        "C4\nLate-Increasing",
    ]
    x = np.arange(len(labels))
    width = 0.34
    ax3.bar(
        x - width / 2,
        align["Python_modal_n"],
        width=width,
        color=colors["Python EM"],
        label="Python EM",
    )
    ax3.bar(
        x + width / 2,
        align["Mplus_modal_n"],
        width=width,
        color=colors["Mplus 7"],
        label="Mplus 7",
    )
    for i, row in align.iterrows():
        y = max(row["Python_modal_n"], row["Mplus_modal_n"]) + 22
        diff = int(row["difference_Mplus_minus_Python"])
        ax3.text(i, y, f"diff {diff:+d}", ha="center", fontsize=8, color="0.25")
    ax3.set_title("C. K=4 modal class counts after alignment", loc="left", fontweight="bold")
    ax3.set_ylabel("Modal-assigned n")
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels)
    ax3.set_ylim(0, max(align["Python_modal_n"].max(), align["Mplus_modal_n"].max()) + 120)
    ax3.grid(axis="y", alpha=0.25)
    ax3.legend(frameon=False, loc="upper left")

    fig.suptitle(
        "Supplementary Figure S1. Mplus replication of the 9-item LCGA solution",
        fontsize=12,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.02,
        "Mplus 7 and Python EM fit indices match to rounding precision; "
        "label-aligned modal agreement = 99.96%, Cohen's kappa = .9994.",
        ha="center",
        fontsize=8,
    )

    for ext, kwargs in {
        "png": {"dpi": 600},
        "tiff": {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}},
        "pdf": {},
    }.items():
        out = OUT / f"FigureS1_Mplus_replication_9item.{ext}"
        fig.savefig(out, bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)
    print(f"Saved Mplus replication figure to {OUT}")


if __name__ == "__main__":
    main()
