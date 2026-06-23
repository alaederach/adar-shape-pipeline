#!/usr/bin/env python3
"""
Rank candidate A->G edit positions from a ShapeMapper profile, to help choose the
read-sorting cutoff for step 02.

Reads a ShapeMapper *_profile.txt, keeps positions whose reference base is the
edit-target base (default A), and writes them sorted by Modified_rate (descending).
Inspect the result (especially for the ADAR-DMSO sample) to pick either a single
min_rate threshold or an explicit list of candidate positions for step 02.
"""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--profile", required=True, help="ShapeMapper *_profile.txt")
    ap.add_argument("--out",     required=True, help="output ranked TSV")
    ap.add_argument("--refbase", default="A",   help="edit-target reference base (default: A)")
    args = ap.parse_args()

    with open(args.profile) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        pos_i  = idx["Nucleotide"]
        base_i = idx["Sequence"]
        rate_i = idx["Modified_rate"]
        depth_i = idx.get("Modified_effective_depth")
        rows = []
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) <= rate_i or c[base_i] != args.refbase:
                continue
            try:
                rate = float(c[rate_i])
            except ValueError:
                continue
            if rate != rate:        # skip nan (positions with no coverage / outside the amplicon)
                continue
            depth = c[depth_i] if depth_i is not None and depth_i < len(c) else ""
            rows.append((int(c[pos_i]), c[base_i], rate, depth))

    rows.sort(key=lambda r: r[2], reverse=True)
    with open(args.out, "w") as o:
        o.write("position\tbase\tmodified_rate\tmodified_effective_depth\n")
        for pos, base, rate, depth in rows:
            o.write("%d\t%s\t%.6f\t%s\n" % (pos, base, rate, depth))

    print("ranked %d %s positions -> %s" % (len(rows), args.refbase, args.out))


if __name__ == "__main__":
    main()
