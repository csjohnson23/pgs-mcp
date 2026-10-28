"""The only numeric cutoffs in the QC layer. Every warning that uses one reports it in evidence."""

# SMALL_SCORE: below this many variants, each variant moves the score noticeably.
SMALL_SCORE_MAX_VARIANTS = 20

# BUILD_MISMATCH_SUSPECTED: median variant match rate across scores below this (percent).
BUILD_MISMATCH_MEDIAN_MATCHED_PCT = 20.0

# Presentation, not QC: full per-score detail is returned for at most this many scores
# unless specific pgs_ids are requested, so tool results stay small enough for a model.
DETAIL_SCORE_LIMIT = 25
