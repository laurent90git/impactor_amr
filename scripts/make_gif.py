"""Animated GIF of a 2D run (schlieren of the density, or another field), one frame per output.

    python scripts/make_gif.py out/gif/impact2009_A -o docs/_static/fig/impact2009_A.gif [--field schlieren|rho|alpha|vonmises]

The time of each frame is read from the run log (<prefix>.log) when the outputs are evenly spaced (--nout).
"""

import argparse
import glob
import io
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from schlieren import COL, raster  # noqa: E402


def frame_image(d, field, k, vlim):
    if field == "schlieren":
        img, ext = raster(d, COL["rho"])
        gy, gx = np.gradient(np.log(img))
        g = np.hypot(gx, gy)
        return np.exp(-k * g / np.nanmax(g)), ext, "gray", (0, 1)
    img, ext = raster(d, COL[field])
    return img, ext, "viridis" if field != "vonmises" else "inferno", vlim


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix", help="output prefix, e.g. out/gif/impact2009_A")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--field", default="schlieren")
    ap.add_argument("--k", type=float, default=40.0)
    ap.add_argument("--tf", type=float, default=None, help="final time (frame k at k tf / nframes)")
    ap.add_argument("--fps", type=float, default=8.0)
    ap.add_argument("--title", default="")
    a = ap.parse_args()

    files = sorted(f for f in glob.glob(a.prefix + "_[0-9][0-9][0-9].npy"))
    if not files:
        sys.exit("no frames for " + a.prefix)
    # times: the run writes frame 0 at t = 0, then nout frames evenly spaced up to tf
    tf = a.tf
    if tf is None:
        log = a.prefix + ".log"
        if os.path.exists(log):
            last = [l for l in open(log) if not l.startswith("#")][-1].split()
            tf = float(last[1])
    nfr = len(files) - 1
    vlim = None
    if a.field != "schlieren":
        lo, hi = np.inf, -np.inf
        for f in files[:: max(1, len(files) // 8)]:
            v = np.load(f)[:, COL[a.field]]
            lo, hi = min(lo, np.nanpercentile(v, 1)), max(hi, np.nanpercentile(v, 99))
        vlim = (lo, hi)
    frames = []
    for i, f in enumerate(files):
        d = np.load(f)
        img, ext, cmap, lim = frame_image(d, a.field, a.k, vlim)
        fig, ax = plt.subplots(figsize=(5.2, 5.2 * (ext[3] - ext[2]) / (ext[1] - ext[0]) + 0.5))
        ax.imshow(img, origin="lower", extent=ext, cmap=cmap, vmin=lim[0] if lim else None, vmax=lim[1] if lim else None)
        t = tf * i / nfr if tf else None
        ax.set_title(f"{a.title}  t = {t * 1e6:7.1f} us" if t is not None else os.path.basename(f), fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=80)
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=128))
    frames[0].save(a.out, save_all=True, append_images=frames[1:], duration=int(1000 / a.fps), loop=0, optimize=True)
    print(f"wrote {a.out} ({len(frames)} frames, {os.path.getsize(a.out) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
