"""Animated GIF of a 2D run (schlieren of the density, or another field), one frame per output.

    python scripts/make_gif.py out/gif/impact2009_A -o docs/_static/fig/impact2009_A.gif [--field schlieren|rho|alpha|vonmises|level|p|u|v]
    python scripts/make_gif.py out/gif/impact2009_A -o a.gif --field schlieren,level     # side-by-side panels

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
from schlieren import COL, apply_slice, raster  # noqa: E402


def frame_image(d, field, k, vlim):
    if field == "schlieren":
        img, ext = raster(d, COL["rho"])
        gy, gx = np.gradient(np.log(img))
        g = np.hypot(gx, gy)
        return np.exp(-k * g / np.nanmax(g)), ext, "gray", (0, 1)
    img, ext = raster(d, COL[field])
    cmap = {"vonmises": "inferno", "level": "tab10"}.get(field, "viridis")
    return img, ext, cmap, vlim


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix", help="output prefix, e.g. out/gif/impact2009_A")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--field", default="schlieren")
    ap.add_argument("--k", type=float, default=40.0)
    ap.add_argument("--tf", type=float, default=None, help="final time (frame k at k tf / nframes)")
    ap.add_argument("--fps", type=float, default=8.0)
    ap.add_argument("--title", default="")
    ap.add_argument("--slice", default="", help="3D: plane to show, e.g. z=0.5")
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
    fields = a.field.split(",")
    vlims = {}
    for fld in fields:
        if fld == "schlieren":
            continue
        lo, hi = np.inf, -np.inf
        for f in files[:: max(1, len(files) // 8)]:
            v = apply_slice(np.load(f), a.slice)[:, COL[fld]]
            lo, hi = min(lo, np.nanpercentile(v, 1)), max(hi, np.nanpercentile(v, 99))
        vlims[fld] = (lo, hi) if fld != "level" else (np.floor(lo), np.ceil(hi))
    frames = []
    for i, f in enumerate(files):
        d = apply_slice(np.load(f), a.slice)
        imgs = [frame_image(d, fld, a.k, vlims.get(fld)) for fld in fields]
        ext = imgs[0][1]
        aspect = (ext[3] - ext[2]) / (ext[1] - ext[0])
        fig, axs = plt.subplots(1, len(fields), figsize=(4.6 * len(fields), 4.6 * aspect + 0.6), squeeze=False)
        for ax, fld, (img, ext, cmap, lim) in zip(axs[0], fields, imgs):
            im = ax.imshow(img, origin="lower", extent=ext, cmap=cmap, vmin=lim[0] if lim else None, vmax=lim[1] if lim else None)
            if fld != "schlieren":
                fig.colorbar(im, ax=ax, shrink=0.75)
            ax.set_xlabel(fld, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
        t = tf * i / nfr if tf else None
        fig.suptitle(f"{a.title}  t = {t * 1e6:7.1f} us  ({len(d)} cells)" if t is not None else os.path.basename(f), fontsize=10)
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
