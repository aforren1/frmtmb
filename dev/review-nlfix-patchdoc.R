f <- "dev/reviews/2026-11-02-nlfix.md"
s <- readChar(f, file.info(f)$size, useBytes = TRUE)
rp <- function(a, b) { if (!grepl(a, s, fixed = TRUE)) stop("missing: ", a); s <<- sub(a, b, s, fixed = TRUE) }
rp("- All L1 on my own specs, 27 of 27 readable, none failing on a body.",
   "- L1 on my own specs (`dev/review-nlfix-specs.R`, 40 with a record): 33 identities. The rest: `r_student_sig` and its linear control `r_student_ctl`, not body related (Filed); `r_compois`, `r_gev`, `r_mix_perloc`, translator gaps with L4 agreeing; `r_powzero`, Minor 11; `r_stan_funs`, my own name clash, replaced by `r_stan_funs2`, an identity. `r_hurdle_disc` has no record: it is the segfault.")
rp("`pow()` is new and documented with \"Stan's values\". The tape gives NaN where Stan gives\n    0.", "`pow()` is new and documented with \"Stan's values\". The tape's gradient is NaN where\n    Stan's is finite.")
rp("- Large response, default start: `a + b * x`, y about 1e6, logLik -3047.75 against -1658.76\n  linear, conv = 0, with \"a step along them raises the log-likelihood\". The same on\n  rellib-r16 (`dev/review-nlfix-a13.R`).",
   "- Large response, default start: `a + b * x`, y about 1e6, logLik -3045.62 against -1658.76\n  linear, conv = 0, with \"a step along them raises the log-likelihood\". The same on\n  rellib-r16 (`dev/review-nlfix-a10.R`, seed 9901).")
rp("| brms posterior, default inits (`sampling`, 2 chains, 500 iter) | | b2 mean 0.74, sd 0.38 |", "| brms posterior, default inits (`sampling`, 2 chains, 500 iter, short: Rhat 1.09) | | b2 mean 0.74, sd 0.38 |")
rp("| **d_cloglog, d_lse, d_lde, d_hypot** | **survive every tier** |", "| **d_cloglog, d_lse, d_lde** | **survive test-nl-family.R and test-nl-stan-funs.R, gated and ungated** |\n| **d_hypot** | **survives both files ungated** (gated not run) |")
rp("Every body with a linear equivalent\nreaches it (gap at most 1e-9)", "Every body with a linear equivalent\nreaches it (gap at most 4e-9)")
rp("frmtmb.sample's 74/79 is the same.", "frmtmb.sample's 74/79 is the same: its 5 are `test-brms-suite-*.R` too.")
rp("-267.12 and -255.52 true.", "-267.12 and -255.52 true (`dev/nlfix-d1check.R` rerun: -267.1211, -255.5164).")
writeChar(s, f, eos = NULL, useBytes = TRUE)
cat("patched\n")
