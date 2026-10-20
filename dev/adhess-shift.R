# Lane adhess, finding 8 of the illse2 review (dev/illse2-rev-misc.R in
# lane illse2's tree, part 1): ML on sleepstudy with Days shifted by s,
# against lme4 (ML) on Days. A shift is an exact reparameterization, so
# the maximum log-likelihood is lme4's and the slope's SE is lme4's.
#   ADHESS_ARM=lane|base Rscript dev/adhess-shift.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
data("sleepstudy", package = "lme4")
set.seed(1)
ssu <- sleepstudy[sort(sample(nrow(sleepstudy), 120)), ]
for (sh in c(0, 1e3, 1e4, 1e5, 1e6)) {
  d <- ssu
  d$t <- d$Days + sh
  r <- cap(frm(Reaction ~ t + (1 | Subject), data = d))
  f <- r$v
  l <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = d, REML = FALSE)
  se_l <- sqrt(diag(as.matrix(vcov(l))))[["Days"]]
  fe <- suppressWarnings(fixef(f))
  se_f <- tryCatch(fe["t", "Est.Error"],
                   error = function(e) NA)
  cat(sprintf(paste0("ML shift %-6g: code %d | logLik - lme4 %+.3e | ",
                     "slope - lme4 %+.2e | SE(slope) / lme4 - 1 %+.2e | ",
                     "lost %s | warn %s\n"),
              sh, f$opt$convergence,
              as.numeric(logLik(f)) - as.numeric(logLik(l)),
              fe["t", "Estimate"] - lme4::fixef(l)[["Days"]],
              se_f / se_l - 1, lost_of(f), wtag(r$w)))
}
