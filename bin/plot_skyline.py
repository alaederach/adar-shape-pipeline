#!/usr/bin/env python3
"""Two-panel skyline control figure for the ADAR-SHAPE deconvolution.

Top panel — reactivity:
    bulk standard-SHAPE (blue) vs the sorted/deconvolved SHAPE-ADAR edited-molecule
    reactivity (orange). Both on the SAME normalization scale (pipeline feeds the
    matched --normout / --scaleout outputs of one normalize call).

Bottom panel — editing (the basis for read sorting):
    ADAR-DMSO mutation rate (red, = A->G editing at A positions) vs the non-ADAR
    DMSO control (grey, background). A horizontal line marks the sorting cutoff and
    vertical lines mark every position the reads were sorted on.

Drawing matches RNAvigate's plot_profile_skyline (ax.plot drawstyle="steps-mid" +
fill_between step="mid"); matplotlib only, RNAvigate is not a dependency.
"""
import argparse
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SHAPE_COLOR     = "#1f77b4"  # bulk SHAPE (top, blue)
ADAR_COLOR      = "#ff7f0e"  # edited SHAPE-ADAR (top, orange)
DMSO_ADAR_COLOR = "#d62728"  # ADAR-DMSO editing (bottom, red)
DMSO_CTRL_COLOR = "#7f7f7f"  # DMSO control background (bottom, grey)
CUTOFF_COLOR    = "#2ca02c"  # cutoff line + sort positions (green)


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def read_cols(path, cols):
    """Read Nucleotide, Sequence and the requested value columns from a profile.txt."""
    nt, seq = [], []
    data = {c: [] for c in cols}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for c in cols:
            if c not in reader.fieldnames:
                raise SystemExit(f"{path}: column '{c}' not found. "
                                 f"Available: {reader.fieldnames}")
        for row in reader:
            nt.append(int(row["Nucleotide"]))
            seq.append((row.get("Sequence") or "").upper())
            for c in cols:
                data[c].append(_f(row.get(c)))
    return (np.array(nt), np.array(seq, dtype=object),
            {c: np.array(v, float) for c, v in data.items()})


def skyline(ax, nt, val, color, label, err=None):
    ax.plot(nt, val, drawstyle="steps-mid", color=color, lw=1.0, label=label)
    if err is not None and np.isfinite(err).any():
        ax.fill_between(nt, val - err, val + err,
                        step="mid", color=color, alpha=0.25, lw=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    # top panel (reactivity)
    ap.add_argument("--shape", required=True, help="bulk standard-SHAPE normalized profile (blue)")
    ap.add_argument("--adar", required=True, help="edited SHAPE-ADAR normalized profile (orange)")
    ap.add_argument("--column", default="Norm_profile")
    ap.add_argument("--errcol", default="Norm_stderr", help="stderr column ('' to disable)")
    ap.add_argument("--shape-label", default="SHAPE (bulk)")
    ap.add_argument("--adar-label", default="SHAPE-ADAR (edited)")
    # bottom panel (editing)
    ap.add_argument("--adar-dmso", required=True, help="ADAR-DMSO profile (editing rate)")
    ap.add_argument("--dmso-control", required=True, help="profile holding the non-ADAR DMSO control")
    ap.add_argument("--rate-col", default="Modified_rate", help="rate column in --adar-dmso")
    ap.add_argument("--dmso-control-col", default="Untreated_rate",
                    help="rate column to read from --dmso-control")
    ap.add_argument("--adar-dmso-label", default="ADAR-DMSO (editing)")
    ap.add_argument("--dmso-control-label", default="DMSO control")
    ap.add_argument("--rate-ylabel", default="Mutation rate")
    ap.add_argument("--cutoff", type=float, default=None,
                    help="sorting min-rate: draws the horizontal line and (if no "
                         "--sort-positions) derives sorted positions as refbase "
                         "positions at/above it")
    ap.add_argument("--sort-positions", default=None,
                    help="explicit comma-separated sorted positions (overrides --cutoff)")
    ap.add_argument("--refbase", default="A", help="reference base sorted on")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    errcol = args.errcol or None
    top_cols = [args.column] + ([errcol] if errcol else [])
    nt_s, _, ds = read_cols(args.shape, top_cols)
    nt_a, _, da = read_cols(args.adar, top_cols)
    nt_ad, seq_ad, dad = read_cols(args.adar_dmso, [args.rate_col])
    nt_dc, _, ddc = read_cols(args.dmso_control, [args.dmso_control_col])
    nmax = max(int(nt_s.max()), int(nt_a.max()), int(nt_ad.max()), int(nt_dc.max()))

    # Positions the reads were sorted on.
    if args.sort_positions:
        sortpos = [int(x) for x in args.sort_positions.replace(" ", "").split(",") if x]
    elif args.cutoff is not None:
        rb = args.refbase.upper().replace("T", "U")
        rate = dad[args.rate_col]
        sortpos = [int(p) for p, s, r in zip(nt_ad, seq_ad, rate)
                   if s.replace("T", "U") == rb and np.isfinite(r) and r >= args.cutoff]
    else:
        sortpos = []

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    # ── top: reactivity ──
    for p in sortpos:  # sorted positions carried up from the editing panel
        ax1.axvline(p, color=CUTOFF_COLOR, lw=0.8, ls=":", alpha=0.6, zorder=0)
    skyline(ax1, nt_s, ds[args.column], SHAPE_COLOR, args.shape_label,
            ds.get(errcol) if errcol else None)
    skyline(ax1, nt_a, da[args.column], ADAR_COLOR, args.adar_label,
            da.get(errcol) if errcol else None)
    ax1.axhline(0, color="0.7", lw=0.6, zorder=0)
    ax1.set_ylabel(args.column.replace("_", " "))
    ax1.legend(frameon=False, loc="upper right", fontsize=8)
    if args.title:
        ax1.set_title(args.title)

    # ── bottom: editing + sorting cutoff/positions ──
    for p in sortpos:  # behind the data
        ax2.axvline(p, color=CUTOFF_COLOR, lw=0.8, ls=":", alpha=0.6, zorder=0)
    if sortpos:  # one legend proxy for the vertical lines
        ax2.plot([], [], color=CUTOFF_COLOR, lw=0.8, ls=":",
                 label=f"sorted positions (n={len(sortpos)})")
    skyline(ax2, nt_ad, dad[args.rate_col], DMSO_ADAR_COLOR, args.adar_dmso_label)
    skyline(ax2, nt_dc, ddc[args.dmso_control_col], DMSO_CTRL_COLOR, args.dmso_control_label)
    if args.cutoff is not None:
        ax2.axhline(args.cutoff, color=CUTOFF_COLOR, lw=1.0, ls="--",
                    label=f"cutoff = {args.cutoff:g}")
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel(args.rate_ylabel)
    ax2.set_xlabel("Nucleotide position")
    ax2.set_xlim(1, nmax)
    ax2.legend(frameon=False, loc="upper right", fontsize=8)

    fig.set_size_inches(max(10, nmax / 25.0), 6.4)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"wrote {args.out}  ({nmax} nt, {len(sortpos)} sorted positions)")


if __name__ == "__main__":
    main()
