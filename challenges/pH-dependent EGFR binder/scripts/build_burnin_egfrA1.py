#!/usr/bin/env python3
"""Burn-in training table for the egfrA1 single-site optimizer, from the egfrW10 campaign.

100 sequences x 10 draws = 1,000 already-folded, already-scored structures -> one row per SEQUENCE.

LABEL:
    label = mean n_epitope over that sequence's GATE-PASSING draws (sc_rmsd < 1.5 A)
    floor = gate_frac >= 0.5, else the sequence is not a valid CANDIDATE
⛔ The floor is a candidate filter, never a term summed into the label (the dftfit lesson).
⚠️ Below-floor sequences are KEPT in the table with passes_floor=0 so the training step can decide
whether to learn from them; dropping data silently is how a burn-in gets quietly cut down.

usage: build_burnin_egfrA1.py <w100.tsv> <gated.csv> <out.csv>
"""
import csv, collections, statistics, sys

GATE = 1.5
FLOOR = 0.5


def main(tsv, gated, outp):
    seqs = {}
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        seqs[f[0].rsplit("_id", 1)[1]] = f[max(range(len(f)), key=lambda i: len(f[i]))]
    per = collections.defaultdict(list)
    for r in csv.DictReader(open(gated)):
        per[r["seq"]].append(dict(nep=int(r["n_epitope"]), sc=float(r["sc_rmsd"]),
                                  tgt=int(r["n_target_contacts"]), gate=int(r["passes_gate"])))
    assert set(per) == set(seqs), (f"⛔ sequence/label mismatch: {len(seqs)} seqs, {len(per)} scored; "
                                   f"symmetric difference {set(per) ^ set(seqs)}")
    L = {len(s) for s in seqs.values()}
    assert len(L) == 1, f"⛔ ragged sequence lengths {L}"
    print(f"{len(seqs)} sequences of {L.pop()} aa, {sum(len(v) for v in per.values())} draws")

    rows = []
    for sid, v in per.items():
        g = [d for d in v if d["gate"]]
        gf = len(g) / len(v)
        rows.append(dict(
            seq_id=f"id{sid}", sequence=seqs[sid], n_draws=len(v), n_gated=len(g),
            gate_frac=round(gf, 3),
            label=round(statistics.mean(d["nep"] for d in g), 4) if g else "",
            mean_nep_all=round(statistics.mean(d["nep"] for d in v), 4),
            full_epitope_frac=round(sum(1 for d in v if d["nep"] == 11) / len(v), 3),
            median_sc=round(statistics.median(d["sc"] for d in v), 3),
            mean_tgt=round(statistics.mean(d["tgt"] for d in v), 1),
            passes_floor=int(gf >= FLOOR)))
    rows.sort(key=lambda r: -(r["label"] if r["label"] != "" else -1))
    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    ok = [r for r in rows if r["passes_floor"]]
    bad = [r for r in rows if not r["passes_floor"]]
    print(f"pass floor (gate_frac >= {FLOOR}): {len(ok)}   below floor: {len(bad)}")
    print(f"label: min {min(r['label'] for r in ok):.2f}  median "
          f"{statistics.median(r['label'] for r in ok):.2f}  max {max(r['label'] for r in ok):.2f}")
    print(f"\nTOP 10 (the set whose consensus is the starting sequence):")
    for r in ok[:10]:
        print(f"  {r['seq_id']:<6} label {r['label']:5.2f}  gate {r['n_gated']}/10  "
              f"full-epi {100*r['full_epitope_frac']:3.0f}%  med_sc {r['median_sc']:.2f}")
    print(f"\nBELOW-FLOOR sequences (kept, flagged):")
    for r in bad:
        lab = r["label"] if r["label"] != "" else "n/a"
        print(f"  {r['seq_id']:<6} gate {r['n_gated']}/10  label_on_gated {lab}  "
              f"mean_all {r['mean_nep_all']:.2f}  med_sc {r['median_sc']:.2f}")
    print(f"\nwrote {outp}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
