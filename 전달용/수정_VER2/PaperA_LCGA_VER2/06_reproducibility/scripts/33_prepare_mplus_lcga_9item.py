# -*- coding: utf-8 -*-
"""
Prepare Mplus replication inputs for the 9-item primary LCGA.

Mplus 7-era executables can be brittle with non-ASCII paths, so the actual
Mplus run directory is an ASCII-only folder:

    D:/mplus_kcyps_paperA

The generated model mirrors the Python/R primary LCGA:
  - y1-y7 = 9-item depression, E04 excluded
  - keep participants with >=3 valid waves
  - quadratic growth factors with time scores -3,-2,-1,0,1,2,3
  - LCGA: growth-factor variances fixed to zero
  - one class-invariant residual variance shared by all seven waves
  - K=1..5 class solutions
"""
from __future__ import annotations

from pathlib import Path
import shutil

import numpy as np
import pandas as pd


ROOT = Path("D:/2026/SCI/청소년패널조사")
DESC = ROOT / "output" / "desc"
P9 = ROOT / "output" / "paperA_lcga" / "primary_9item"
MPLUS_DIR = Path("D:/mplus_kcyps_paperA")


def clean_mplus_dir() -> None:
    MPLUS_DIR.mkdir(parents=True, exist_ok=True)
    for pattern in ("lcga9_k*.inp", "lcga9_k*.out", "lcga9_k*_cprob.dat", "*.gh5"):
        for p in MPLUS_DIR.glob(pattern):
            p.unlink()


def write_data() -> pd.DataFrame:
    long_df = pd.read_parquet(DESC / "long_panel.parquet")
    y = long_df.pivot_table(index="ID", columns="wave", values="depression9")[list(range(1, 8))]
    valid_count = y.notna().sum(axis=1)
    y = y.loc[valid_count >= 3].copy()
    y.columns = [f"y{w}" for w in range(1, 8)]

    out = y.reset_index()
    out["ID"] = out["ID"].astype(int)
    cols = ["ID"] + [f"y{w}" for w in range(1, 8)]
    out = out[cols].replace({np.nan: -999})

    out.to_csv(
        MPLUS_DIR / "kcyps_dep9.dat",
        sep=" ",
        header=False,
        index=False,
        float_format="%.6f",
    )
    out.to_csv(P9 / "mplus_input_dep9_wide.csv", index=False, encoding="utf-8-sig")
    return out


def model_block(k: int) -> str:
    lines = [
        "MODEL:",
        "  %OVERALL%",
        "    i s q | y1@-3 y2@-2 y3@-1 y4@0 y5@1 y6@2 y7@3;",
        "    i@0;",
        "    s@0;",
        "    q@0;",
        "    y1-y7 (resid);",
    ]
    if k > 1:
        for c in range(1, k + 1):
            lines.extend([
                f"  %c#{c}%",
                "    [i s q];",
            ])
    return "\n".join(lines)


def write_inp(k: int) -> None:
    savedata = ""
    if k == 4:
        savedata = """
SAVEDATA:
  FILE IS lcga9_k4_cprob.dat;
  SAVE = CPROBABILITIES;
"""

    text = f"""TITLE:
  KCYPS 2018 m1 9-item depression LCGA K={k};

DATA:
  FILE IS kcyps_dep9.dat;

VARIABLE:
  NAMES ARE id y1 y2 y3 y4 y5 y6 y7;
  USEVARIABLES ARE y1 y2 y3 y4 y5 y6 y7;
  IDVARIABLE IS id;
  MISSING ARE ALL (-999);
  CLASSES = c({k});

ANALYSIS:
  TYPE = MIXTURE;
  ESTIMATOR = ML;
  STARTS = 1000 250;
  STITERATIONS = 20;
  PROCESSORS = 4;

{model_block(k)}

OUTPUT:
  TECH1 TECH4 TECH8;
{savedata}
"""
    (MPLUS_DIR / f"lcga9_k{k}.inp").write_text(text, encoding="ascii")


def main() -> None:
    clean_mplus_dir()
    out = write_data()
    for k in range(1, 6):
        write_inp(k)

    readme = f"""Mplus 9-item LCGA replication inputs

Run directory: {MPLUS_DIR.as_posix()}
Data: kcyps_dep9.dat
N analytic sample: {len(out)}
Variables: id y1-y7
Missing: -999
Model: quadratic LCGA, class-invariant shared residual variance, i/s/q variances fixed to zero.

Suggested run command:
  & 'C:\\Program Files\\Mplus\\Mplus.exe' lcga9_k1.inp
  ...
  & 'C:\\Program Files\\Mplus\\Mplus.exe' lcga9_k5.inp
"""
    (MPLUS_DIR / "README.txt").write_text(readme, encoding="ascii")

    project_dir = P9 / "mplus"
    project_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MPLUS_DIR / "README.txt", project_dir / "README.txt")
    print(f"Wrote Mplus data and K=1..5 inputs to {MPLUS_DIR}")
    print(f"Analytic N = {len(out)}")


if __name__ == "__main__":
    main()
