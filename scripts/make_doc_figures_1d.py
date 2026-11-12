"""1D verification runs and figures of the documentation (docs/_static/fig). Runs the C++ solver (one core).

    python scripts/make_doc_figures_1d.py [--skip-runs]
"""

import argparse
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "prototype"))
import hvi1d  # noqa: E402

EXE = os.path.join(ROOT, "build", "hvi1d")
OUT = os.path.join(ROOT, "out", "doc1d")
FIG = os.path.join(ROOT, "docs", "_static", "fig")
C = {"x": 0, "h": 3, "rho": 5, "alpha": 6, "u": 7, "v": 8, "p": 10, "s11": 11, "s12": 13, "a21": 14}


def run(case, tag, *args):
    cmd = [EXE, "--case", case, "--out", OUT, "--tag", tag, "--log-every", "100000000", *map(str, args)]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def load(case, tag, k=0):
    a = np.load(os.path.join(OUT, f"{case}_{tag}_{k:03d}.npy"))
    return a[np.argsort(a[:, 0])]


def panel(axs, a, keys, style="k.-", label=None, xscale=1.0):
    for ax, k in zip(axs, keys):
        ax.plot(a[:, 0] * xscale, a[:, C[k]], style, ms=2, lw=0.7, label=label)
        ax.set_title(k)
        ax.grid(alpha=0.3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-runs", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    cases = ["advect", "slip-advect", "shear-slab", "double-expansion", "solid-shear-tube", "solid-gas-tube", "impact1d"]
    keys = ["rho", "alpha", "u", "v", "p", "s11", "s12", "a21"]
    n = 640  # level 5 on [0, 1]

    if not a.skip_runs:
        for c in cases:
            run(c, "l5", "--level", 5)
        run("slip-advect", "noslip", "--level", 5, "--alim", 0)
        run("shear-slab", "l7", "--level", 7)
        run("shear-slab", "late", "--level", 7, "--tf", 180e-6)
        run("impact1d", "uni8", "--level", 8)
        run("impact1d", "mr8", "--level", 8, "--coarse-level", 3)
        run("kanel1d", "l9", "--level", 9)
        run("spall1d", "mr", "--level", 10, "--coarse-level", 4, "--paper-times", "--nout", 5)

    # ---- ladder vs prototype -----------------------------------------------------------------------------------------
    rows = []
    for c in cases:
        r = load(c, "l5")
        model = hvi1d.Model([hvi1d.COPPER, hvi1d.AIR], alim=1e-2)
        x, al, rl, v, p, tf, _ = hvi1d.setup(c, n)
        U = model.run(hvi1d.initial(model, x, al, rl, v, p), 1.0 / n, tf, 0.4, 2, verbose=False)
        f = hvi1d.fields(model, U)
        d = max(np.abs(r[:, C[k]] - f[k]).mean() / max(np.abs(f[k]).max(), 1e-30) for k in ("rho", "u", "v", "p"))
        rows.append((c, tf, d))
        fig, axs = plt.subplots(2, 4, figsize=(15, 6))
        for ax, k in zip(axs.ravel(), keys):
            ax.plot(x, f[k], "r-", lw=1.2, label="NumPy prototype")
        panel(axs.ravel(), r, keys, "k.", "C++ (samurai)")
        axs[0, 0].legend(fontsize=8)
        fig.suptitle(f"{c}, 640 cells, t = {tf:.3e} s")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, f"ladder_{c}.png"), dpi=90)
        plt.close(fig)
    with open(os.path.join(FIG, "ladder_table.rst"), "w") as fo:
        fo.write(".. list-table::\n   :header-rows: 1\n\n   * - case\n     - final time [s]\n     - max rel. L1 difference (rho, u, v, p)\n")
        for c, tf, d in rows:
            fo.write(f"   * - ``{c}``\n     - {tf:.3e}\n     - {d:.1e}\n")

    # ---- slip: with and without fluidisation --------------------------------------------------------------------------
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6))
    for tag, st, lab in (("l5", "k.-", "fluidisation, alim = 1e-2"), ("noslip", "r.-", "no slip treatment")):
        panel(axs, load("slip-advect", tag), ["v", "p", "s12", "a21"], st, lab)
    axs[0].legend(fontsize=8)
    fig.suptitle("Favrie 2009 sec. 6.1: interface advected at 1000 m/s with a tangential velocity jump of 100 m/s")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "slip_advect.png"), dpi=90)
    plt.close(fig)

    # ---- shear slab vs the linear solution ----------------------------------------------------------------------------
    ct = np.sqrt(hvi1d.COPPER.mu / hvi1d.COPPER.rho0)
    fig, axs = plt.subplots(2, 3, figsize=(14, 7))
    for row, (tag, t) in enumerate((("l7", 54.6e-6), ("late", 180e-6))):
        r = load("shear-slab", tag)
        xx = np.linspace(0, 1, 4000)
        vex, aex = hvi1d.dalembert_slab(xx, t, ct)
        panel(axs[row], r, ["v", "a21", "p"])
        axs[row, 0].plot(xx, vex, "r-", lw=1, label="linear elasticity")
        axs[row, 1].plot(xx, aex, "r-", lw=1)
        axs[row, 0].set_ylabel(f"t = {t * 1e6:.1f} us")
        axs[row, 0].legend(fontsize=8)
    fig.suptitle("Favrie 2009 sec. 6.2: shear waves in a copper slab with free surfaces, 2560 cells")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "shear_slab.png"), dpi=90)
    plt.close(fig)

    # ---- multiresolution vs uniform -----------------------------------------------------------------------------------
    uni, mr = load("impact1d", "uni8"), load("impact1d", "mr8")
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6))
    panel(axs[:3], uni, ["rho", "u", "p"], "r-", f"uniform level 8 ({len(uni)} cells)")
    panel(axs[:3], mr, ["rho", "u", "p"], "k.", f"multiresolution levels 3-8 ({len(mr)} cells)")
    axs[3].plot(mr[:, 0], np.log2(1.0 / (mr[:, 3] / 0.05)), "k.", ms=2)
    axs[3].set_title("level")
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "mr_impact1d.png"), dpi=90)
    plt.close(fig)

    # ---- Kanel / Ndanou: elastic precursors and plastic plateau --------------------------------------------------------
    k = load("kanel1d", "l9")
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6))
    panel(axs, k, ["rho", "u", "p", "s11"], "k.-", None, 1e3)
    for ax in axs:
        ax.set_xlim(1.5, 4.6)
        ax.set_xlabel("x [mm]")
    fig.suptitle("Ndanou 2015 sec. 6.2.1: Al (700 m/s) onto Ti, t = 0.2 us, 8704 cells")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "kanel1d.png"), dpi=90)
    plt.close(fig)

    def at(xmm, col):
        return k[np.argmin(np.abs(k[:, 0] * 1e3 - xmm)), C[col]]

    with open(os.path.join(FIG, "kanel_table.rst"), "w") as fo:
        fo.write(".. list-table::\n   :header-rows: 1\n\n   * - quantity\n     - analytical (Ndanou 2015)\n     - this code\n")
        vals = [
            ("Al elastic precursor: density", "2715.51", f"{at(1.70, 'rho'):.2f}"),
            ("Al elastic precursor: velocity", "691.64", f"{at(1.70, 'u'):.2f}"),
            ("Ti elastic precursor: density", "4581.90", f"{at(4.10, 'rho'):.2f}"),
            ("Ti elastic precursor: velocity", "74.58", f"{at(4.10, 'u'):.2f}"),
            ("Ti elastic precursor: pressure [GPa]", "1.41", f"{at(4.10, 'p') / 1e9:.3f}"),
            ("plastic plateau: Al density", "2880.90", f"{at(2.5, 'rho'):.2f}"),
            ("plastic plateau: Ti density", "4782.80", f"{at(3.5, 'rho'):.2f}"),
            ("plastic plateau: interface velocity", "297.149", f"{at(3.0, 'u'):.2f}"),
            ("plastic plateau: normal stress [GPa]", "-7.50285", f"{at(3.0, 's11') / 1e9:.4f}"),
        ]
        for q, an, nu in vals:
            fo.write(f"   * - {q}\n     - {an}\n     - {nu}\n")

    # ---- spallation ---------------------------------------------------------------------------------------------------
    files = sorted(f for f in os.listdir(OUT) if f.startswith("spall1d_mr_") and f.endswith(".npy"))
    if files:
        fig, ax = plt.subplots(figsize=(10, 3.6))
        for f in files:
            s = np.load(os.path.join(OUT, f))
            s = s[np.argsort(s[:, 0])]
            ax.plot(s[:, 0] * 1e3, s[:, 5], lw=0.8, label=f)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("density")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
        fig.suptitle("Ndanou 2015 sec. 6.2.2: spallation of the Ti plate (yield 1030 MPa)")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "spall1d.png"), dpi=90)
        plt.close(fig)
    print("figures in", FIG)


if __name__ == "__main__":
    main()
