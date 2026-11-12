"""Numerical schlieren (and other fields) of 2D outputs, rasterised on the finest level.

    python scripts/schlieren.py out/impact2009_001.npy [more.npy ...] -o fig.png [--field schlieren|rho|p|vonmises|alpha|level]
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

COL = {"rho": 5, "alpha": 6, "u": 7, "v": 8, "p": 10, "vonmises": 16, "level": 4, "chi": 15}


def raster(a, col):
    x, y, h = a[:, 0], a[:, 1], a[:, 3]
    hmin = h.min()
    x0, y0 = (x - 0.5 * h).min(), (y - 0.5 * h).min()
    nx = int(round(((x + 0.5 * h).max() - x0) / hmin))
    ny = int(round(((y + 0.5 * h).max() - y0) / hmin))
    img = np.full((ny, nx), np.nan)
    i0 = np.round((x - 0.5 * h - x0) / hmin).astype(int)
    j0 = np.round((y - 0.5 * h - y0) / hmin).astype(int)
    k = np.round(h / hmin).astype(int)
    for kk in np.unique(k):
        m = k == kk
        for di in range(kk):
            for dj in range(kk):
                img[j0[m] + dj, i0[m] + di] = a[m, col]
    return img, (x0, x0 + nx * hmin, y0, y0 + ny * hmin)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--out", default="schlieren.png")
    ap.add_argument("--field", default="schlieren")
    ap.add_argument("--k", type=float, default=40.0)
    a = ap.parse_args()

    n = len(a.files)
    ncol = min(n, 4)
    nrow = (n + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.2 * nrow), squeeze=False)
    for ax, f in zip(axs.ravel(), a.files):
        d = np.load(f)
        if a.field == "schlieren":
            img, ext = raster(d, COL["rho"])
            gy, gx = np.gradient(np.log(img))  # log: the air shocks are visible next to the copper/air jump
            g = np.hypot(gx, gy)
            img = np.exp(-a.k * g / np.nanmax(g))
            ax.imshow(img, origin="lower", extent=ext, cmap="gray", vmin=0, vmax=1)
        else:
            img, ext = raster(d, COL[a.field])
            im = ax.imshow(img, origin="lower", extent=ext, cmap="viridis")
            fig.colorbar(im, ax=ax, shrink=0.8)
        ax.set_title(f"{os.path.basename(f)}  ({len(d)} cells)", fontsize=9)
        ax.set_aspect("equal")
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(a.out, dpi=110)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
