"""Compare a 1D run of the C++ solver with the NumPy prototype on the same mesh.

    python scripts/compare_proto.py <case> <run.npy> [--n 320] [--order 2] [--alim 1e-2] [--plot out.png]

Prints the relative L1 difference of rho, u, v, p, s11, s12 (normalised by the max of the prototype field).
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "prototype"))
import hvi1d  # noqa: E402

COLS = {"rho": 5, "alpha_s": 6, "u": 7, "v": 8, "p": 10, "s11": 11, "s12": 13, "a21": 14}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("npy")
    ap.add_argument("--n", type=int, default=320)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--alim", type=float, default=1e-2)
    ap.add_argument("--cfl", type=float, default=0.4)
    ap.add_argument("--tf", type=float, default=None)
    ap.add_argument("--plot", default=None)
    a = ap.parse_args()

    run = np.load(a.npy)
    run = run[np.argsort(run[:, 0])]
    model = hvi1d.Model([hvi1d.COPPER, hvi1d.AIR], alim=a.alim, slip=a.alim > 0)
    x, al, rl, v, p, tf, _ = hvi1d.setup(a.case, a.n)
    if a.tf is not None:
        tf = a.tf
    U = hvi1d.initial(model, x, al, rl, v, p)
    U = model.run(U, 1.0 / a.n, tf, a.cfl, a.order, verbose=False)
    f = hvi1d.fields(model, U)
    if run.shape[0] != a.n:
        # adapted run: compare on the prototype's mesh by sampling the cell containing each point
        idx = np.clip(np.searchsorted(run[:, 0] - 0.5 * run[:, 3], x, side="right") - 1, 0, run.shape[0] - 1)
        run = run[idx]
    worst = 0.0
    for k, c in COLS.items():
        ref = f[k]
        scale = max(np.abs(ref).max(), 1e-30)
        d = np.abs(run[:, c] - ref).mean() / scale
        worst = max(worst, d) if k in ("rho", "u", "v", "p") else worst
        print(f"  {k:8s} rel L1 diff {d:.3e}   (C++ range [{run[:, c].min():.4e}, {run[:, c].max():.4e}], proto [{ref.min():.4e}, {ref.max():.4e}])")
    print(f"  worst (rho, u, v, p): {worst:.3e}")
    if a.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2, 4, figsize=(16, 7))
        for ax, (k, c) in zip(axs.ravel(), COLS.items()):
            ax.plot(x, f[k], "r-", lw=1.0, label="prototype")
            ax.plot(run[:, 0], run[:, c], "k.", ms=2, label="C++")
            ax.set_title(k)
            ax.grid(alpha=0.3)
        axs[0, 0].legend()
        fig.suptitle(f"{a.case}: C++ vs prototype, n={a.n}")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=100)
    return worst


if __name__ == "__main__":
    main()
