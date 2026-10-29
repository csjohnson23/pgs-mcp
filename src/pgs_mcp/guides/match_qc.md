# Interpreting pgsc_calc match QC

> This guide encodes domain judgement. Contributions from PGS practitioners are welcome.

This guide is for interpreting the output of the `assess_match_qc` tool. The tool reports
raw numbers and structural warnings. Your job is to explain what the numbers mean for *this*
score in *this* target dataset, and how much they should change a user's confidence in the
calculated scores.

## Principles

1. **Report, then interpret.** Quote the numbers you rely on. Separate the data from what you
   infer about cause or impact.
2. **No universal cutoffs.** A given match rate means different things for a 20-variant
   clinical score and a 1M-variant genome-wide score. Reason from score size, weight
   distribution and the pattern of missingness, rather than from a single threshold.
3. **The pipeline verdict is one input.** `pipeline_pass` only says the variant match rate
   cleared `min_overlap`. It says nothing about which variants were lost.
4. **Distinguish run-level from score-level insights.** If more than one score exists in the
   log, a problem shared by every score (a missing chromosome, no indels) is a property of the
   target data, and should be reported as such.
5. **Limit interpretation to known parameters.** This QC describes coverage, not score
   validity. A perfectly matched score can still be poorly constructed or poorly suited to the
   target population, but the values here cannot ascertain that.
6. **The target cohort is unknown to the tool.** Ancestry, genotyping platform, imputation
   panel and quality filters are not in pgsc_calc's matching outputs. If the user has told you
   about them, use that context and say so; if an interpretation depends on them, ask rather
   than assume.

## Inputs

Check `log_used` first. **Without the match log**, only variant-level shares from the summary
are available: weight coverage, the largest unmatched variants, per-chromosome coverage, indel
matching and missing positions are all unknown (`LOG_NOT_SUPPLIED`). Say so, recommend a
re-run with the log included, and limit interpretation: do not describe a score as well
reproduced on variant match rate alone, since a score can match 70% of its variants and keep
only 16% of its weight.

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
  still removes real between-person variance.
- **Weight coverage uses |effect weight| alone, without allele frequency.** A variant's
  contribution to score variance scales with weight² × 2p(1−p), so a lost rare variant with a
  large weight costs more weight coverage than variance, and a lost common variant of moderate
  weight the reverse. Treat weight coverage as a proxy for retained signal, not a measure of
  it, and say so when a conclusion rests on it.

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

A/T and C/G variants cannot be strand-resolved by allele comparison. By default pgsc_calc
excludes them. With `keep_ambiguous` on it keeps them, and because direct matches always
outrank flipped ones, every ambiguous variant is scored in the unflipped orientation: it
assumes the target and scoring file share a strand. That is usually safe for well-harmonised
data on the forward strand, and risky for raw array data. A mis-stranded ambiguous variant
contributes its effect with the wrong sign, silently. Weigh the share of ambiguous variants
and whether the target data is known to be strand-aligned (ask if it matters and is unknown).

### Strand flips (`flipped_pct`)

Variants that matched only after complementing alleles. Ambiguous variants are never matched
by a flip, so this rate describes non-ambiguous variants only. Moderate rates are normal for
array data. Very high rates in imputed or sequenced data suggest a strand convention problem.

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
| `MISSING_POSITION` | Variants with no chromosome/position | Entries with no genomic coordinate: HLA alleles, presence/absence codes (e.g. P/N), copy-number or other non-SNV entries | These cannot be matched by position; judge impact from weight coverage and `largest_unmatched`, and note them as unscorable in this pipeline |
| `DUPLICATES_EXCLUDED` | Variants dropped for duplicate IDs | Duplicate records in the target | Deduplicate target variants |
| `AMBIGUOUS_RETAINED` | Ambiguous variants were kept | `keep_ambiguous` on | See strand-ambiguous section |
| `SMALL_SCORE` | Few variants in the score | Score design | Every missing variant matters; read the variant list directly |
| `BUILD_MISMATCH_SUSPECTED` | Uniformly low matching | Genome build or ID convention mismatch | Confirm target build matches the scoring files |
| `LOG_SUMMARY_MISMATCH` | Log and summary disagree | Files from different runs | Re-run QC with files from one run before interpreting |
| `LOG_NOT_SUPPLIED` | Only the summary was given | User omitted the log | See "What was supplied"; weight-based conclusions are impossible |

## Writing the interpretation

**For one score:** lead with a one-sentence bottom line on how faithfully the score was
reproduced in this dataset. Then cover, only where relevant: weight vs variant coverage,
what was lost and why (with the warning codes), strand handling, contextualized warnings, and
concrete next steps.

**For many scores:** lead with run-level findings that affect every score and their likely
upstream cause. Then group scores by pattern (e.g. "passed but lost most of their weight",
"failed due to indel-heavy design") rather than walking through scores one by one. Name the
scores most in need of attention and why. Note where the inferred settings differ from
pgsc_calc defaults (`min_overlap` 0.75, ambiguous and multiallelic variants excluded).

Keep it to what a user needs to decide whether, and how, to use these scores.
