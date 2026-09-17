## The names variables(fit) and the draws columns carry, checked against
## the names brms ITSELF would give the same model, derived from brms's
## own parse rather than from frmtmb's labeler. dev/brmsnames-match.R
## cannot check this: its shim is built from frmtmb's names, so a
## variables() row there compares a name with itself.
##
##   Rscript dev/brmsnames-naming.R base > dev/brmsnames-log/naming-base.txt
##   Rscript dev/brmsnames-naming.R lane > dev/brmsnames-log/naming-lane.txt
##
## brms's expected names come from brms::default_prior() (one row per
## population-level coefficient and group-level SD, with dpar, nlpar and
## resp) and brms's own brms:::combine_prefix() for the prefix, and from
## the empty brmsfit's ranef frame for the r_ names. Data seed 41.
arm <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/brmsnames-libs.R")
brmsnames_libs(arm)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
q(requireNamespace("brms"))
cat("arm", arm, " frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n\n")

set.seed(41)
n <- 240
dd <- data.frame(x = stats::rnorm(n), z = stats::runif(n),
                 g = factor(rep(1:12, 20)))
dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + stats::rnorm(12, 0, 0.5)[dd$g],
                     exp(0.2 * dd$z))
dd$y1 <- dd$y
dd$y2 <- stats::rnorm(n, dd$x, 1)
dd$yn <- 2 * exp(-0.7 * dd$z) + stats::rnorm(n, 0, 0.1)
# hostile names (punch round 1, MAJOR 1): a factor level with a space
# and one with a hyphen, a group whose levels carry whitespace, an
# interaction group, and responses with an underscore and a dot
dd$f <- factor(rep(c("a b", "c-d", "e"), length.out = n))
dd$gs <- factor(paste("lev", rep(1:8, length.out = n)))
dd$h <- factor(rep(c("p", "q", "r"), length.out = n))
dd$y_a <- dd$y
dd$y.b <- dd$y2

brms_expected <- function(bform, data, family) {
  pr <- q(brms::default_prior(bform, data = data, family = family))
  pr <- as.data.frame(pr)
  pre <- function(r) {
    brms:::combine_prefix(list(resp = r$resp, dpar = r$dpar,
                               nlpar = r$nlpar))
  }
  usc <- function(p, s) if (nzchar(p)) paste0(p, "_", s) else s
  out <- character(0)
  # a smooth's unpenalized columns are listed under class b, and brms's
  # rename_sm() names them bs_<prefix>_<label>_<k>
  sm <- pr[pr$class == "sds" & nzchar(pr$coef), , drop = FALSE]
  sm_lab <- brms:::rename(sm$coef)
  for (i in seq_len(nrow(pr))) {
    r <- pr[i, ]
    is_sm <- r$class == "b" && nzchar(r$coef) &&
      any(startsWith(r$coef, paste0(sm_lab, "_")))
    if (is_sm) {
      out <- c(out, paste0("bs_", usc(pre(r), r$coef)))
    } else if (r$class == "sds" && nzchar(r$coef)) {
      out <- c(out, paste0("sds_", usc(pre(r), brms:::rename(r$coef)), "_1"))
    } else if (r$class %in% c("sigma", "shape", "nu", "phi", "zi", "hu") &&
               !nzchar(r$coef) && !nzchar(r$dpar) && !nzchar(r$nlpar)) {
      # a distributional parameter nobody wrote a formula for is the
      # parameter itself, sigma or sigma_<resp>
      out <- c(out, if (nzchar(r$resp)) paste0(r$class, "_", r$resp) else
        r$class)
    } else if (r$class == "b" && nzchar(r$coef)) {
      out <- c(out, paste0("b_", usc(pre(r), r$coef)))
    } else if (r$class == "Intercept" && !nzchar(r$coef)) {
      out <- c(out, paste0("b_", usc(pre(r), "Intercept")))
    } else if (r$class == "sd" && nzchar(r$coef) && nzchar(r$group)) {
      out <- c(out, paste0("sd_", r$group, "__", usc(pre(r), r$coef)))
    }
  }
  unique(out)
}

cases <- list(
  list("y ~ x + (1 + x | g)",
       frmtmb::bf(y ~ x + (1 + x | g)), brms::bf(y ~ x + (1 + x | g)),
       gaussian()),
  list("bf(y ~ x + (1 | g), sigma ~ z + (1 | g))",
       frmtmb::bf(y ~ x + (1 | g), sigma ~ z + (1 | g)),
       brms::bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), gaussian()),
  list("mvbf(y1 ~ x + (1 | g), y2 ~ x + (1 | g), sigma ~ z), no rescor",
       frmtmb::mvbf(frmtmb::bf(y1 ~ x + (1 | g), sigma ~ z),
                    frmtmb::bf(y2 ~ x + (1 | g), sigma ~ z), rescor = FALSE),
       brms::mvbf(brms::bf(y1 ~ x + (1 | g), sigma ~ z),
                  brms::bf(y2 ~ x + (1 | g), sigma ~ z), rescor = FALSE),
       gaussian()),
  list("yn ~ a * exp(-b * z), a ~ 1 + (1 | g), b ~ 1, nl",
       frmtmb::bf(yn ~ a * exp(-b * z), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
       brms::bf(yn ~ a * exp(-b * z), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
       gaussian()),
  list("y ~ x + I(x^2) + poly(z, 2) + f + (1 | gs) + (1 | g:h), hostile",
       frmtmb::bf(y ~ x + I(x^2) + poly(z, 2) + f + (1 | gs) + (1 | g:h)),
       brms::bf(y ~ x + I(x^2) + poly(z, 2) + f + (1 | gs) + (1 | g:h)),
       gaussian()),
  list("mvbf(y_a ~ x + (1 | g), y.b ~ f), responses with _ and .",
       frmtmb::mvbf(frmtmb::bf(y_a ~ x + (1 | g)), frmtmb::bf(y.b ~ f),
                    rescor = FALSE),
       brms::mvbf(brms::bf(y_a ~ x + (1 | g)), brms::bf(y.b ~ f),
                  rescor = FALSE),
       gaussian()),
  list("y ~ s(x) + (1 | g), smooth",
       frmtmb::bf(y ~ s(x) + (1 | g)), brms::bf(y ~ s(x) + (1 | g)),
       gaussian()),
  list("negbinomial yc ~ x, shape unmodeled",
       frmtmb::bf(yc ~ x), brms::bf(yc ~ x), negbinomial(),
       brms::negbinomial()),
  # every other unmodeled dpar the review asked about (MAJOR 3)
  list("student y ~ x + (1 | g), sigma and nu unmodeled",
       frmtmb::bf(y ~ x + (1 | g)), brms::bf(y ~ x + (1 | g)), student(),
       brms::student()),
  list("Beta yb ~ x, phi unmodeled",
       frmtmb::bf(yb ~ x), brms::bf(yb ~ x), Beta(), brms::Beta()),
  list("zero_inflated_poisson yc ~ x, zi unmodeled",
       frmtmb::bf(yc ~ x), brms::bf(yc ~ x), zero_inflated_poisson(),
       brms::zero_inflated_poisson()),
  list("hurdle_gamma yh ~ x, shape and hu unmodeled",
       frmtmb::bf(yh ~ x), brms::bf(yh ~ x), hurdle_gamma(),
       brms::hurdle_gamma()),
  list("student bf(y ~ x, nu ~ 1), nu written out",
       frmtmb::bf(y ~ x, nu ~ 1), brms::bf(y ~ x, nu ~ 1), student(),
       brms::student())
)
dd$yb <- stats::plogis(0.3 * dd$x + stats::rnorm(n, 0, 0.5))
dd$yh <- ifelse(stats::runif(n) < 0.3, 0, stats::rgamma(n, 2, 2))
dd$yc <- stats::rpois(n, exp(1 + 0.3 * dd$x))

tot_expected <- 0L
tot_found <- 0L
for (cs in cases) {
  cat("== ", cs[[1L]], " ==\n", sep = "")
  fam_b <- if (length(cs) >= 5L) cs[[5L]] else cs[[4L]]
  want <- tryCatch(brms_expected(cs[[3L]], dd, fam_b),
                   error = function(e) paste("ERROR", conditionMessage(e)))
  fam <- cs[[4L]]
  fit <- tryCatch(q(frm(cs[[2L]], family = fam, data = dd,
                        start = if (grepl("nl", cs[[1L]]))
                          list(beta = c(1, 1)))),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat("frm() failed:", conditionMessage(fit), "\n\n")
    next
  }
  got <- variables(fit)
  hit <- want %in% got
  tot_expected <- tot_expected + length(want)
  tot_found <- tot_found + sum(hit)
  cat("brms names:          ", paste(want, collapse = " "), "\n")
  cat("variables(fit):      ", paste(got, collapse = " "), "\n")
  cat("brms names missing:  ", paste(want[!hit], collapse = " "), "\n")
  # default_prior() has no row for a correlation, so cor_ and rescor__
  # names are not in brms's list; anything else extra is a name brms
  # does not have
  extra <- setdiff(got[!grepl("^(cor_|rescor__)", got)], want)
  cat("names brms lacks:    ", paste(extra, collapse = " "), "\n")
  cat("found ", sum(hit), " of ", length(want), "\n", sep = "")
  # the r_ names of the draws columns, against brms's ranef frame
  # the base build has no brms_par_labels(); its draws are labeled by
  # frmtmb.sample's own labeler, which is what its draws columns carry
  labs <- tryCatch(frmtmb::brms_par_labels(fit), error = function(e) {
    tryCatch(frmtmb.sample:::all_par_labels(fit),
             error = function(e2) NULL)
  })
  b <- q(brms::brm(cs[[3L]], data = dd, family = fam_b, empty = TRUE,
                   backend = "rstan"))
  re <- b$ranef
  if (!is.null(labs) && NROW(re)) {
    rw <- character(0)
    for (i in seq_len(nrow(re))) {
      p <- brms:::combine_prefix(re[i, ])
      # brms's own levels (combine_groups() joins an interaction with
      # _), with rename_re_levels()'s whitespace-to-dot
      lev <- gsub("[[:space:]]", ".", attr(re, "levels")[[re$group[i]]])
      rw <- c(rw, paste0("r_", re$group[i],
                         if (nzchar(p)) paste0("__", p), "[", lev, ",",
                         re$coef[i], "]"))
    }
    rhit <- rw %in% labs
    tot_expected <- tot_expected + length(rw)
    tot_found <- tot_found + sum(rhit)
    cat("r_ names found in the draws labels: ", sum(rhit), " of ",
        length(rw), "\n", sep = "")
    if (any(!rhit)) {
      cat("  missing, first 3: ",
          paste(utils::head(rw[!rhit], 3L), collapse = " "), "\n")
    }
  } else if (is.null(labs)) {
    cat("r_ names: this build has no brms_par_labels()\n")
  }
  cat("\n")
}
cat("== total, arm ", arm, ": brms names present ", tot_found, " of ",
    tot_expected, " ==\n", sep = "")
cat("DONE\n")
