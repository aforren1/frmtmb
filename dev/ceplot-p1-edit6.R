# Lane ceplot punch 1 (m5): posterior_samples(pars = ) orders the
# population coefficients by brms's rule, not by fixef()'s.
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/extensions/frmtmb.sample/R/methods-draws.R"
s <- readLines(f)
a <- grep("^#' The variables `posterior_samples[(]pars = [)]` selects, in brms's order[.]$", s)
b <- grep("^ps_select <- function", s)
e <- b + which(s[(b + 1):length(s)] == "}")[1]
stopifnot(length(a) == 1, length(b) == 1)
new <- c(
"#' The variables `posterior_samples(pars = )` selects, in brms's order.",
"#'",
"#' brms's `extract_pars()`: each pattern keeps the variables it matches",
"#' in `variables()` order, and `fixed = TRUE` keeps the named ones in",
"#' the order given. brms's `variables()` lists the intercept of every",
"#' distributional parameter first (`b_Intercept`, `b_sigma_Intercept`,",
"#' an ordinal fit's `b_Intercept[k]`), then every other population",
"#' coefficient in the order of its predictor, a nonlinear parameter's",
"#' intercept staying with its own coefficients (`b_a_Intercept`,",
"#' `b_a_z`, `b_b_Intercept`), as measured on brms 2.23.0",
"#' (`dev/ceplot-p1-psorder-brms.R`). frmtmb's `variables()` lists each",
"#' predictor's coefficients together, so the `b_` columns are put in",
"#' brms's order before the patterns are matched. `variables()` itself",
"#' keeps its own order, so that `posterior_samples(x)` with no `pars`",
"#' still has the columns of `variables(x)`, as brms's does.",
"#'",
"#' @noRd",
"ps_select <- function(x, pars, fixed) {",
"  v <- variables(x)",
"  bi <- grep(\"^b_\", v)",
"  if (length(bi) > 1L) {",
"    fit <- draws_base_fit(x)",
"    nlp <- unique(unlist(lapply(fit$spec$responses, `[[`, \"nlpars\")))",
"    resp <- names(fit$spec$responses)",
"    is_nl <- if (length(nlp)) {",
"      grepl(paste0(\"^b_((\", paste(resp, collapse = \"|\"), \")_)?(\",",
"                   paste(nlp, collapse = \"|\"), \")_\"), v[bi])",
"    } else {",
"      rep(FALSE, length(bi))",
"    }",
"    front <- grepl(\"Intercept([[][0-9]+[]])?$\", v[bi]) & !is_nl",
"    v[bi] <- c(v[bi][front], v[bi][!front])",
"  }",
"  draws_extract_pars(pars, v, fixed)",
"}")
s <- c(s[seq_len(a - 1)], new, s[(e + 1):length(s)])
writeLines(s, f)
cat("edited\n")
