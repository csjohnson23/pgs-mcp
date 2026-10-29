# Match QC: EXAMPLE

- **Scores:** 3
- **Pipeline verdict:** 3 pass, 0 fail
- **Genome build:** GRCh38
- **Variant log:** used

## Run-wide warnings

- **`CHROM_ABSENT`**: No score variant found in the target on chromosome(s) X (10) (score-variant pairs in parentheses). The target data likely lacks these chromosomes or labels them differently.
- **`TARGET_LACKS_INDELS`**: 0 of 7 indel/non-SNV score variants matched while 153 SNVs did. The target data likely contains no indels (array data, or indels filtered or not normalised).

## Inferred pipeline settings

Inferred from the match summary, not read from the run configuration:

- `min_overlap`: at most 0.800 (no score failed)
- `keep_ambiguous`: yes
- `keep_multiallelic`: unknown

## Distribution across scores

| Metric | Min | Median | Max |
|---|---|---|---|
| Variants matched (%) | 80.00 | 82.00 | 98.33 |
| Weight coverage (%) | 25.29 | 47.75 | 98.76 |
| Weight minus variant coverage (points) | -54.71 | -34.25 | 0.43 |
| Variants per score | 15 | 60 | 100 |

## Chromosome coverage (score-variant pairs)

| Chromosome | Variants | Matched | Matched % |
|---|---|---|---|
| 1 | 8 | 8 | 100.0% |
| 2 | 8 | 7 | 87.5% |
| 3 | 12 | 11 | 91.7% |
| 4 | 6 | 6 | 100.0% |
| 5 | 11 | 10 | 90.9% |
| 6 | 9 | 9 | 100.0% |
| 7 | 2 | 2 | 100.0% |
| 8 | 9 | 7 | 77.8% |
| 9 | 4 | 4 | 100.0% |
| 10 | 4 | 4 | 100.0% |
| 11 | 8 | 7 | 87.5% |
| 12 | 9 | 9 | 100.0% |
| 13 | 13 | 12 | 92.3% |
| 14 | 3 | 3 | 100.0% |
| 15 | 8 | 8 | 100.0% |
| 16 | 13 | 12 | 92.3% |
| 17 | 9 | 8 | 88.9% |
| 18 | 4 | 4 | 100.0% |
| 19 | 3 | 3 | 100.0% |
| 20 | 6 | 6 | 100.0% |
| 21 | 4 | 4 | 100.0% |
| 22 | 9 | 9 | 100.0% |
| X | 10 | 0 | 0.0% |

## Warnings by score

| Warning | Scores |
|---|---|
| `INDELS_UNMATCHED` | 2 |
| `AMBIGUOUS_RETAINED` | 1 |
| `CHROM_ABSENT` | 1 |
| `SMALL_SCORE` | 1 |
| `MISSING_POSITION` | 1 |

## All scores

Sorted by weight coverage, lowest first.

| Score | Variants | Pipeline | Matched | Weight coverage | Warnings |
|---|---|---|---|---|---|
| hla_small_hmPOS_GRCh38 | 15 | pass | 80.0% | 25.3% | SMALL_SCORE, INDELS_UNMATCHED, MISSING_POSITION |
| large_effect_loss_hmPOS_GRCh38 | 100 | pass | 82.0% | 47.8% | AMBIGUOUS_RETAINED, CHROM_ABSENT, INDELS_UNMATCHED |
| mostly_clean_hmPOS_GRCh38 | 60 | pass | 98.3% | 98.8% |  |

## Score details

### large_effect_loss_hmPOS_GRCh38

#### Metrics

| Metric | Value |
|---|---|
| Variants in score | 100 |
| Matched | 82.0% |
| Matched, unambiguous | 76.0% |
| Strand-ambiguous | 6.0% |
| Matched after strand flip | 4.0% |
| Excluded | 0.0% |
| Unmatched | 18.0% |
| Weight coverage (share of total abs. weight matched) | 47.8% |
| Largest abs. weight in score | 0.4791 |
| Largest unmatched (effect/other, abs. weight) | 2:3437810 TA/T, 0.4791; 16:15264840 A/AG, 0.4253; 11:70019441 CAT/C, 0.396 |
| SNVs matched | 82 / 96 |
| Indels/non-SNVs matched | 0 / 4 |
| Variants missing a position | 0 |

#### Warnings

- **`AMBIGUOUS_RETAINED`**: 6.0% of variants are strand-ambiguous (A/T or C/G) and were kept; their effect allele cannot be confirmed by strand.
- **`CHROM_ABSENT`**: No variants found in the target on chromosome(s) X (10) (score variants in parentheses).
- **`INDELS_UNMATCHED`**: None of the score's 4 indel/non-SNV variants matched.

#### Observations

- Weight coverage is 47.8% against a variant match rate of 82.0% (-34.2 points): the variants that did not match carry more than their share of total |effect weight|.
- The largest-effect variant in the score (2:3437810 TA/T, |weight| 0.4791, 9.8% of total |weight|) is unmatched.
- Variant match rate per autosome ranges from 66.7% (chr13) to 100.0% (chr9) across 22 autosomes.
- 4 of 18 unmatched or excluded variants are indels/non-SNVs.

#### Match types

| Match type | Variants |
|---|---|
| refalt | 78 |
| refalt_flip | 4 |
| discarded candidate matches (not_best) | 0 |

#### Chromosome coverage

| Chromosome | Variants | Matched | Matched % |
|---|---|---|---|
| 1 | 5 | 5 | 100.0% |
| 2 | 7 | 6 | 85.7% |
| 3 | 7 | 6 | 85.7% |
| 4 | 2 | 2 | 100.0% |
| 5 | 8 | 7 | 87.5% |
| 6 | 6 | 6 | 100.0% |
| 7 | 1 | 1 | 100.0% |
| 8 | 4 | 3 | 75.0% |
| 9 | 1 | 1 | 100.0% |
| 10 | 2 | 2 | 100.0% |
| 11 | 6 | 5 | 83.3% |
| 12 | 4 | 4 | 100.0% |
| 13 | 3 | 2 | 66.7% |
| 14 | 2 | 2 | 100.0% |
| 15 | 4 | 4 | 100.0% |
| 16 | 9 | 8 | 88.9% |
| 17 | 6 | 5 | 83.3% |
| 18 | 3 | 3 | 100.0% |
| 19 | 2 | 2 | 100.0% |
| 20 | 2 | 2 | 100.0% |
| 21 | 2 | 2 | 100.0% |
| 22 | 4 | 4 | 100.0% |
| X | 10 | 0 | 0.0% |

### mostly_clean_hmPOS_GRCh38

#### Metrics

| Metric | Value |
|---|---|
| Variants in score | 60 |
| Matched | 98.3% |
| Matched, unambiguous | 98.3% |
| Strand-ambiguous | 0.0% |
| Matched after strand flip | 0.0% |
| Excluded | 0.0% |
| Unmatched | 1.7% |
| Weight coverage (share of total abs. weight matched) | 98.8% |
| Largest abs. weight in score | 0.03 |
| Largest unmatched (effect/other, abs. weight) | 8:98966210 A/C, 0.0152 |
| SNVs matched | 59 / 60 |
| Indels/non-SNVs matched | 0 / 0 |
| Variants missing a position | 0 |

#### Warnings

No score-level warnings.

#### Observations

- Weight coverage is 98.8% against a variant match rate of 98.3% (+0.4 points): the variants that did not match carry less than their share of total |effect weight|.
- Variant match rate per autosome ranges from 75.0% (chr8) to 100.0% (chr9) across 21 autosomes.

#### Match types

| Match type | Variants |
|---|---|
| refalt | 59 |
| discarded candidate matches (not_best) | 0 |

#### Chromosome coverage

| Chromosome | Variants | Matched | Matched % |
|---|---|---|---|
| 1 | 1 | 1 | 100.0% |
| 2 | 1 | 1 | 100.0% |
| 3 | 3 | 3 | 100.0% |
| 4 | 4 | 4 | 100.0% |
| 6 | 3 | 3 | 100.0% |
| 7 | 1 | 1 | 100.0% |
| 8 | 4 | 3 | 75.0% |
| 9 | 2 | 2 | 100.0% |
| 10 | 2 | 2 | 100.0% |
| 11 | 2 | 2 | 100.0% |
| 12 | 5 | 5 | 100.0% |
| 13 | 9 | 9 | 100.0% |
| 14 | 1 | 1 | 100.0% |
| 15 | 3 | 3 | 100.0% |
| 16 | 4 | 4 | 100.0% |
| 17 | 3 | 3 | 100.0% |
| 18 | 1 | 1 | 100.0% |
| 19 | 1 | 1 | 100.0% |
| 20 | 4 | 4 | 100.0% |
| 21 | 1 | 1 | 100.0% |
| 22 | 5 | 5 | 100.0% |

### hla_small_hmPOS_GRCh38

#### Metrics

| Metric | Value |
|---|---|
| Variants in score | 15 |
| Matched | 80.0% |
| Matched, unambiguous | 80.0% |
| Strand-ambiguous | 0.0% |
| Matched after strand flip | 0.0% |
| Excluded | 0.0% |
| Unmatched | 20.0% |
| Weight coverage (share of total abs. weight matched) | 25.3% |
| Largest abs. weight in score | 0.9 |
| Largest unmatched (effect/other, abs. weight) | no position, effect HLA-DRB1*15:01, 0.9; no position, effect HLA-DQB1*06:02, 0.7; no position, effect HLA-B*08:01, 0.4 |
| SNVs matched | 12 / 12 |
| Indels/non-SNVs matched | 0 / 3 |
| Variants missing a position | 3 |

#### Warnings

- **`SMALL_SCORE`**: Score has only 15 variants; each one carries substantial weight, so single missing variants shift scores noticeably.
- **`INDELS_UNMATCHED`**: None of the score's 3 indel/non-SNV variants matched.
- **`MISSING_POSITION`**: Variants missing a chromosome or position in the scoring file, which cannot be matched by position: 3.

#### Observations

- Weight coverage is 25.3% against a variant match rate of 80.0% (-54.7 points): the variants that did not match carry more than their share of total |effect weight|.
- The largest-effect variant in the score (no position, effect HLA-DRB1*15:01, |weight| 0.9, 34% of total |weight|) is unmatched.
- All 3 unmatched or excluded variants are indels/non-SNVs.

#### Match types

| Match type | Variants |
|---|---|
| refalt | 12 |
| discarded candidate matches (not_best) | 0 |

#### Chromosome coverage

| Chromosome | Variants | Matched | Matched % |
|---|---|---|---|
| 1 | 2 | 2 | 100.0% |
| 3 | 2 | 2 | 100.0% |
| 5 | 3 | 3 | 100.0% |
| 8 | 1 | 1 | 100.0% |
| 9 | 1 | 1 | 100.0% |
| 13 | 1 | 1 | 100.0% |
| 15 | 1 | 1 | 100.0% |
| 21 | 1 | 1 | 100.0% |

