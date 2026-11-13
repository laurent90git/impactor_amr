"""2D figures of the documentation from the production runs (out/prod) -> docs/_static/fig.

    python scripts/make_doc_figures_2d.py
"""

import glob
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PROD = os.path.join(ROOT, "out", "prod")
FIG = os.path.join(ROOT, "docs", "_static", "fig")
PY = sys.executable
SCH = os.path.join(ROOT, "scripts", "schlieren.py")


def fig(pattern, name, field="schlieren"):
    files = sorted(glob.glob(os.path.join(PROD, pattern)))
    if not files:
        print("no files for", pattern)
        return
    subprocess.run([PY, SCH, *files, "-o", os.path.join(FIG, name), "--field", field], check=True)


def log_summary(name):
    f = os.path.join(PROD, name + ".log")
    if not os.path.exists(f):
        return None
    a = np.loadtxt(f)
    return a


def main():
    os.makedirs(FIG, exist_ok=True)
    for tag in ("A", "B", "C"):
        fig(f"impact2009_{tag}_0*.npy", f"impact2009_{tag}.png")
        fig(f"impact2009_{tag}_0*.npy", f"impact2009_{tag}_rho.png", "rho")
    fig("impact2009_A_0*.npy", "impact2009_A_level.png", "level")
    fig("impact2012_plastic_0*.npy", "impact2012.png")
    fig("impact2012_plastic_0*.npy", "impact2012_rho.png", "rho")
    fig("impact2012_plastic_0*.npy", "impact2012_alpha.png", "alpha")
    fig("penetration2015_0*.npy", "penetration2015_rho.png", "rho")
    fig("penetration2015_0*.npy", "penetration2015_alpha.png", "alpha")

    # conservation table from the logs
    rows = []
    for name in ("impact2009_A", "impact2009_B", "impact2009_C", "impact2012_plastic", "penetration2015"):
        a = log_summary(name)
        if a is None or a.ndim < 2 or len(a) < 2:
            continue
        nph = a.shape[1] - 10  # columns: step t dt cells masses(N) momentum(3) energy floored yielded
        m0, m1 = a[0, 4 : 4 + nph].sum(), a[-1, 4 : 4 + nph].sum()
        e0, e1 = a[0, 4 + nph + 3], a[-1, 4 + nph + 3]
        rows.append((name, a[-1, 1], int(a[-1, 3]), int(a[:, 3].max()), abs(m1 - m0) / m0, (e1 - e0) / abs(e0), int(a[-1, -2])))
    with open(os.path.join(FIG, "impacts_table.rst"), "w") as fo:
        fo.write(".. list-table::\n   :header-rows: 1\n\n")
        fo.write("   * - run\n     - final time [s]\n     - cells (final)\n     - cells (max)\n     - rel. mass change\n"
                 "     - rel. energy change\n     - floored cell-stages\n")
        for r in rows:
            fo.write(f"   * - ``{r[0]}``\n     - {r[1]:.3e}\n     - {r[2]}\n     - {r[3]}\n     - {r[4]:.1e}\n     - {r[5]:.1e}\n     - {r[6]}\n")
    print("done")


if __name__ == "__main__":
    main()
