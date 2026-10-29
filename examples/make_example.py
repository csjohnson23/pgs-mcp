"""Write the synthetic pgsc_calc match outputs in this folder: `python examples/make_example.py`.

The target mimics imputed array data: no chrX and no indels. Three made-up scores:

large_effect_loss  100 variants, 82% matched (passes min_overlap 0.75). The unmatched variants
                   include every chrX variant and all four indels, which carry the largest weights.
mostly_clean       60 variants, 59 matched. Weight coverage tracks variant coverage.
hla_small          15 variants, 80% matched. Its three HLA allele entries have no position and
                   hold most of the score's weight.
"""

import csv
import random
from pathlib import Path

HERE = Path(__file__).parent
DATASET = "EXAMPLE"
UNAMBIGUOUS = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"), ("A", "C"), ("G", "T")]
AMBIGUOUS = [("A", "T"), ("C", "G")]

LOG_COLUMNS = [
    "row_nr", "accession", "chr_name", "chr_position", "effect_allele", "other_allele",
    "effect_weight", "effect_type", "ID", "REF", "ALT", "matched_effect_allele", "match_type",
    "is_multiallelic", "ambiguous", "match_flipped", "best_match", "exclude",
    "duplicate_best_match", "duplicate_ID", "match_IDs", "match_status", "dataset",
]  # fmt: skip
SUMMARY_COLUMNS = [
    "dataset", "accession", "score_pass", "match_status", "ambiguous", "is_multiallelic",
    "duplicate_best_match", "duplicate_ID", "match_flipped", "match_IDs", "count", "percent",
]  # fmt: skip


def variant(rng, chrom, pos, weight, status, ambiguous=False, flipped=False, alleles=None):
    ea, oa = alleles or rng.choice(AMBIGUOUS if ambiguous else UNAMBIGUOUS)
    return dict(
        chr_name=chrom, chr_position=pos, effect_allele=ea, other_allele=oa,
        effect_weight=round(weight * rng.choice((1, -1)), 4),
        status=status, ambiguous=ambiguous, flipped=flipped,
    )  # fmt: skip


def autosomal(rng, n, scale, status, **kw):
    return [
        variant(rng, str(rng.randint(1, 22)), rng.randint(10**5, 2 * 10**8),
                rng.uniform(0.5, 1.5) * scale, status, **kw)
        for _ in range(n)
    ]  # fmt: skip


def build_scores(rng):
    large = (
        autosomal(rng, 72, 0.03, "matched")
        + autosomal(rng, 6, 0.03, "matched", ambiguous=True)
        + autosomal(rng, 4, 0.03, "matched", flipped=True)
        + autosomal(rng, 4, 0.03, "unmatched")
        + [variant(rng, "X", rng.randint(10**6, 1.5 * 10**8), 0.08 * rng.uniform(0.8, 1.2),
                   "unmatched") for _ in range(10)]
        + [variant(rng, str(c), rng.randint(10**6, 10**8), 0.4 * rng.uniform(0.8, 1.2),
                   "unmatched", alleles=a)
           for c, a in [(2, ("TA", "T")), (8, ("G", "GC")), (11, ("CAT", "C")), (16, ("A", "AG"))]]
    )  # fmt: skip
    clean = autosomal(rng, 59, 0.02, "matched") + autosomal(rng, 1, 0.02, "unmatched")
    hla = autosomal(rng, 12, 0.05, "matched") + [
        variant(rng, "", "", w, "unmatched", alleles=(a, ""))
        for a, w in [("HLA-DRB1*15:01", 0.9), ("HLA-DQB1*06:02", 0.7), ("HLA-B*08:01", 0.4)]
    ]  # fmt: skip
    return {
        "large_effect_loss_hmPOS_GRCh38": large,
        "mostly_clean_hmPOS_GRCh38": clean,
        "hla_small_hmPOS_GRCh38": hla,
    }


def flag(x, matched):
    """Unmatched rows leave the match flags blank, as pgsc_calc does."""
    return str(x).lower() if matched else ""


def log_row(i, accession, v):
    matched = v["status"] == "matched"
    variant_id = f"{v['chr_name']}:{v['chr_position']}:{v['other_allele']}:{v['effect_allele']}"
    return dict(
        row_nr=i, accession=accession, effect_type="additive", dataset=DATASET,
        chr_name=v["chr_name"], chr_position=v["chr_position"],
        effect_allele=v["effect_allele"], other_allele=v["other_allele"],
        effect_weight=v["effect_weight"],
        ID=variant_id if matched else "",
        REF=v["other_allele"] if matched else "", ALT=v["effect_allele"] if matched else "",
        matched_effect_allele=v["effect_allele"] if matched else "",
        match_type=("refalt_flip" if v["flipped"] else "refalt") if matched else "",
        is_multiallelic=flag(False, matched), ambiguous=flag(v["ambiguous"], matched),
        match_flipped=flag(v["flipped"], matched), best_match=flag(True, matched),
        exclude=flag(False, matched),
        duplicate_best_match=flag(False, matched), duplicate_ID=flag(False, matched),
        match_IDs="NA" if matched else "", match_status=v["status"],
    )  # fmt: skip


def summary_rows(accession, variants):
    n = len(variants)
    passed = sum(v["status"] == "matched" for v in variants) / n >= 0.75
    groups = {}
    for v in variants:
        matched = v["status"] == "matched"
        key = ("matched", v["ambiguous"], v["flipped"]) if matched else ("unmatched",)
        groups[key] = groups.get(key, 0) + 1
    for key, count in groups.items():
        matched = key[0] == "matched"
        yield dict(
            dataset=DATASET, accession=accession, score_pass=str(passed).lower(),
            match_status=key[0], ambiguous=flag(matched and key[1], matched),
            is_multiallelic=flag(False, matched), duplicate_best_match=flag(False, matched),
            duplicate_ID=flag(False, matched), match_flipped=flag(matched and key[2], matched),
            match_IDs="NA" if matched else "", count=count, percent=100 * count / n,
        )  # fmt: skip


def main():
    scores = build_scores(random.Random(7))
    with (HERE / "example_log.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, LOG_COLUMNS)
        writer.writeheader()
        i = 0
        for accession, variants in scores.items():
            for v in variants:
                writer.writerow(log_row(i, accession, v))
                i += 1
    with (HERE / "example_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, SUMMARY_COLUMNS)
        writer.writeheader()
        for accession, variants in scores.items():
            writer.writerows(summary_rows(accession, variants))


if __name__ == "__main__":
    main()
