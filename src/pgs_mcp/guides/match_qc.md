# Interpreting pgsc_calc match QC

> **Status: draft.** This guide encodes domain judgement and is meant to be reviewed and
> extended by PGS practitioners. Contributions welcome.

This guide is for interpreting the output of the `assess_match_qc` tool. The tool reports
raw numbers and structural warnings; it never grades a score. Your job is to explain what the
numbers mean for *this* score in *this* target dataset, and how much they should change a
user's confidence in the calculated scores.

## Principles

1. **Report, then interpret.** Quote the numbers you rely on. Separate what the data shows
   from what you infer about its cause.
2. **No universal cutoffs.** A given match rate means different things for a 20-variant
   clinical score and a 1M-variant genome-wide score. Reason from score size, weight
   distribution and the pattern of missingness, not from a single threshold.
3. **The pipeline verdict is one input.** `pipeline_pass` only says the variant match rate
   cleared `min_overlap`. It says nothing about which variants were lost.
4. **Distinguish run-level from score-level causes.** A problem shared by every score (a
   missing chromosome, no indels) is a property of the target data. Fix it once, upstream.
5. **Say what is unknown.** Match QC describes coverage, not validity. A perfectly matched
   score can still be poorly suited to the target population.

## What was supplied

Check `log_used` first. **Without the match log**, only variant-level shares from the summary
are available: weight coverage, the largest unmatched variants, per-chromosome coverage, indel
matching and missing positions are all unknown (`LOG_NOT_SUPPLIED`). Say so plainly, and do
not describe a score as well reproduced on variant match rate alone: a score can match 70% of
its variants and keep 16% of its weight. Recommend re-running with the log for any score the
user intends to rely on.

Whole-run calls list every score in `score_index` (a table: `columns` names each position in
`rows`) and leave `scores` empty for large runs; request specific scores to get their full
detail. The run-wide chromosome table is returned only for whole-run calls.

## Metrics

### Variant match rate (`matched_pct`)

Share of the score's variants found in the target data. This is what pgsc_calc compares to
`min_overlap`. On its own it treats every variant as equally important, which they are not.

### Weight coverage (`weight_coverage_pct`)

Share of the score's total |effect weight| carried by matched variants. This is closer to
"how much of the score was actually computed" than the variant match rate.

- **Weight coverage well below variant match rate** means the missing variants are
  disproportionately large-effect. Common in scores built from a few strong GWAS hits plus
  many small effects, and in scores containing HLA, APOE-region or rare large-effect variants
  that arrays and imputation handle poorly. This is the most important pattern to catch: the
  pipeline may pass the score while most of its signal is gone.
- **Weight coverage above variant match rate** means the lost variants are small-effect.
  Usually benign for genome-wide (LDpred/PRS-CS style) scores.
- Missing weight shifts every individual's score by roughly the same *expected* amount, but
  removes real between-person variance. Expect attenuated discrimination and distorted
  absolute risk if the score is used on its original scale.

### Largest unmatched variants (`largest_unmatched`)

The unmatched variants with the largest |effect weight|, with chromosome, position and
effect/other alleles. When the score's single largest effect is among them, say which variant
it is and how much of the score it carries relative to `max_abs_weight`. If you know the trait,
note whether the position falls in a well-known locus for it (HLA on 6p21, APOE on 19q13, for
example), but present that as a lead for the user to confirm, not a finding. Variants with no
position and multi-character alleles are often HLA alleles, copy-number or other non-SNV
entries that positional matching cannot handle. For small scores a single lost variant can
dominate the result.

### Strand-ambiguous variants (`ambiguous_pct`, `AMBIGUOUS_RETAINED`)

A/T and C/G variants cannot be strand-resolved by allele comparison. With `keep_ambiguous`
on, pgsc_calc keeps them, assuming target and scoring file share a strand. That is usually
safe for well-harmonised, imputed data on the forward strand, and risky for raw array data. A
mis-stranded ambiguous variant contributes its effect with the wrong sign. Weigh the share of
ambiguous variants and whether the target data was strand-aligned.

### Strand flips (`flipped_pct`)

Variants that matched only after complementing alleles. Moderate rates are normal for array
data. Very high rates in imputed or sequenced data suggest a strand convention problem.

### Match types

- `refalt` / `altref`: direct allele match; `*_flip`: matched after a strand flip.
- `no_oa_*`: the scoring file has no other allele, so the match relied on one allele only.
  Less certain, particularly at multiallelic sites.
- `not_best` candidates: alternative matches discarded in favour of a better one. Large
  numbers point to multiallelic sites or split records in the target; when `keep_multiallelic`
  is unknown or off, check whether those sites were handled as intended.
- A large score matched almost entirely by one match type (e.g. all `altref`) is normal when
  the target's REF/ALT orientation is consistent; it is not a problem in itself.

### Even SNV loss across chromosomes

A score can fail with no warnings at all: no indels, no ambiguous variants, no missing
positions, just SNVs lost at a similar rate on every autosome. That pattern points to the
scoring file's variant set being poorly covered by the target (imputation panel, post-imputation
INFO/MAF filters, or a score built on a denser or different panel), not to a build mismatch,
which would produce near-zero matching everywhere. Compare per-chromosome match rates in the
score's detail, or the per-autosome range observation, to confirm the loss is even, and say
the cause is inferred.

### Exclusions (`excluded_by_reason`)

- `duplicate_ID` / `duplicate_best_match`: the same target variant matched several score
  variants, or a score variant had tied best matches. pgsc_calc drops them rather than guess.
- `is_multiallelic`, `ambiguous`: dropped by run settings.

## Warnings

| Code | What it means | Typical cause | What to suggest |
|---|---|---|---|
| `CHROM_ABSENT` (run) | No score variant matched on a chromosome | Target lacks the chromosome (often X/Y/MT), or labels it differently (`23` vs `X`, `chr` prefix). Scoring files may themselves mix `X` and `23` | Check target contigs and labels; impute or include sex chromosomes if the score depends on them |
| `CHROM_ABSENT` (score) | This score lost every variant on a chromosome | Usually the run-level cause; if not run-wide, sparse coverage of a region | Check `absent_run_wide` in the evidence before blaming the score |
| `TARGET_LACKS_INDELS` | No indel matched anywhere | Array data, indels filtered out, or indels not normalised/left-aligned | Normalise indels, or accept and quantify the loss per score |
| `INDELS_UNMATCHED` | This score's indels all failed | Usually `TARGET_LACKS_INDELS` | Look at weight coverage to judge impact |
| `MISSING_POSITION` | Variants with no chromosome/position | Scoring file not harmonised for this build | Use the harmonised (`hmPOS`) file for the target build |
| `DUPLICATES_EXCLUDED` | Variants dropped for duplicate IDs | Duplicate records in the target | Deduplicate target variants |
| `AMBIGUOUS_RETAINED` | Ambiguous variants were kept | `keep_ambiguous` on | See strand-ambiguous section |
| `SMALL_SCORE` | Few variants in the score | Score design | Every missing variant matters; read the variant list directly |
| `BUILD_MISMATCH_SUSPECTED` | Uniformly low matching | Genome build or ID convention mismatch | Confirm target build matches the scoring files |
| `LOG_SUMMARY_MISMATCH` | Log and summary disagree | Files from different runs | Re-run QC with files from one run before interpreting |
| `LOG_NOT_SUPPLIED` | Only the summary was given | User omitted the log | See "What was supplied"; weight-based conclusions are impossible |

## Writing the interpretation

**For one score:** lead with a one-sentence bottom line on how faithfully the score was
reproduced in this dataset. Then cover, only where relevant: weight vs variant coverage,
what was lost and why (with the warning codes), strand handling, and concrete next steps.

**For many scores:** lead with run-level findings that affect every score and their likely
upstream cause. Then group scores by pattern (e.g. "passed but lost most of their weight",
"failed due to indel-heavy design") rather than walking through scores one by one. Name the
scores most in need of attention and why. Note where the inferred settings differ from
pgsc_calc defaults (`min_overlap` 0.75, ambiguous and multiallelic variants excluded).

Keep it to what a user needs to decide whether, and how, to use these scores.
