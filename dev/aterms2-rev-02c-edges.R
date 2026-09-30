# Reviewer, claim 2: subset() edge cases and the post-fit methods of a
# multivariate subset model. Seed 223. Log: dev/aterms2-rev-log-02c.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) suppressWarnings(suppressMessages(expr))
show <- function(label, expr) {
  r <- tryCatch(q(expr), error = function(e) paste("ERROR:",
                                                   conditionMessage(e)))
  cat(sprintf("%-34s ", label))
  if (is.character(r) && length(r) == 1L) cat(substr(r, 1, 230), "\n")
  else {
    cat("class", class(r)[1], " dim/len",
        paste(if (is.null(dim(r))) length(r) else dim(r), collapse = "x"),
        "\n")
  }
  invisible(r)
}
brms_try <- function(label, bform, data) {
  r <- tryCatch(q(brms::standata(bform, data = data)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  brms %-29s ", label))
  if (is.character(r)) cat(substr(r, 1, 200), "\n")
  else cat("N_y1", r$N_y1 %||% r$N, "\n")
}
set.seed(223)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n),
                s = rep(c(TRUE, FALSE), n / 2))
d$y1 <- 1 + d$x + rnorm(n)
d$y2 <- d$z + rnorm(n)
d$sF <- FALSE
d$sT <- TRUE
d$sNA <- d$s
d$sNA[5] <- NA
d$s01 <- as.numeric(d$s)
d$s02 <- ifelse(d$s, 2, 0)
d$sch <- ifelse(d$s, "TRUE", "FALSE")
d$sfac <- factor(d$sch)
d$syes <- ifelse(d$s, "yes", "no")

cat("== subset variable edge cases (frmtmb, then brms standata)\n")
mk <- function(sv) {
  bquote(frm(bf(y1 | subset(.(as.name(sv))) ~ x) + bf(y2 ~ z), data = d,
             family = gaussian()))
}
ref <- q(frm(bf(y1 | subset(s) ~ x) + bf(y2 ~ z), data = d,
             family = gaussian()))
for (sv in c("sF", "sT", "sNA", "s01", "s02", "sch", "sfac", "syes")) {
  f <- show(sv, eval(mk(sv)))
  if (inherits(f, "frmtmb_fit") && sv != "sT") {
    cat("    rows y1", length(f$frame$y$y1), " fn equal to logical s:",
        identical(f$obj$fn(ref$opt$par), ref$obj$fn(ref$opt$par)), "\n")
  }
  if (inherits(f, "frmtmb_fit") && sv == "sT") {
    f0 <- q(frm(bf(y1 ~ x) + bf(y2 ~ z), data = d, family = gaussian()))
    set.seed(1)
    p <- f0$opt$par + rnorm(length(f0$opt$par), 0, 0.1)
    cat("    all TRUE: fn identical to no subset:",
        identical(f$obj$fn(p), f0$obj$fn(p)),
        " gr identical:", identical(f$obj$gr(p), f0$obj$gr(p)), "\n")
    show("    all TRUE fitted() no resp", fitted(f))
  }
  brms_try(sv, brms::bf(as.formula(paste0("y1 | subset(", sv, ") ~ x"))) +
             brms::bf(y2 ~ z) + brms::set_rescor(FALSE), d)
}
show("subset(x > 0) expression", frm(bf(y1 | subset(x > 0) ~ x) +
                                       bf(y2 ~ z), data = d,
                                     family = gaussian()))
brms_try("subset(x > 0)", brms::bf(y1 | subset(x > 0) ~ x) +
           brms::bf(y2 ~ z) + brms::set_rescor(FALSE), d)
show("subset(TRUE) constant", frm(bf(y1 | subset(TRUE) ~ x) + bf(y2 ~ z),
                                  data = d, family = gaussian()))
show("univariate all FALSE", frm(y1 | subset(sF) ~ x, data = d))
show("univariate na.exclude", frm(y1 | subset(s) ~ x, data = d,
                                  na.action = stats::na.exclude))
show("rescor TRUE", frm(bf(y1 | subset(s) ~ x) + bf(y2 ~ z) +
                          set_rescor(TRUE), data = d, family = gaussian()))
show("mv ar()", frm(bf(y1 | subset(s) ~ x + ar(p = 1)) + bf(y2 ~ z),
                    data = d, family = gaussian()))
d$sx <- d$x + rnorm(n, 0, 0.1)
d$sdx <- 0.1
show("me() refused", frm(bf(y1 | subset(s) ~ me(sx, sdx)) + bf(y2 ~ z),
                         data = d, family = gaussian()))

cat("\n== post-fit methods of a multivariate subset model\n")
f <- ref
show("print", capture.output(print(f)))
show("summary", summary(f))
show("nobs", nobs(f))
cat("   nobs value", nobs(f), "\n")
show("AIC", AIC(f))
show("logLik", logLik(f))
show("vcov", vcov(f))
show("confint", confint(f))
show("fixef", fixef(f))
show("ranef", ranef(f))
show("fitted() no resp", fitted(f))
show("fitted(resp = 'y1')", fitted(f, resp = "y1"))
show("fitted(resp = c('y1','y2'))", fitted(f, resp = c("y1", "y2")))
show("predict() no resp", predict(f))
show("predict(resp = 'y1')", predict(f, resp = "y1", ndraws = 10))
show("frm_linpred() no resp", frm_linpred(f))
show("residuals()", residuals(f))
show("residuals(resp = 'y1')", residuals(f, resp = "y1"))
show("simulate()", simulate(f, nsim = 1, seed = 1))
show("conditional_effects()", conditional_effects(f))
show("conditional_effects(resp)", conditional_effects(f, resp = "y1"))
show("hypothesis", hypothesis(f, "y1_x > 0"))
show("model.frame", model.frame(f))
show("update(newdata)", update(f, newdata = d[1:40, ]))
show("vcov_cluster", vcov_cluster(f, cluster = seq_len(n)))
show("dharma_residuals", dharma_residuals(f))
show("pp_check", pp_check(f))
show("emmeans", emmeans::emmeans(f, ~ x, resp = "y1"))
nd <- d[1:10, ]
show("fitted(newdata, resp y1)", fitted(f, newdata = nd, resp = "y1"))
cat("   rows expected", sum(nd$s), "\n")
show("fitted(newdata w/o s)", fitted(f, newdata = nd[, c("x", "z")],
                                     resp = "y1"))
show("fitted(newdata, resp y2)", fitted(f, newdata = nd, resp = "y2"))
cat("\n== univariate subset post-fit\n")
fu <- q(frm(y1 | subset(s) ~ x, data = d))
show("nobs", nobs(fu))
cat("   nobs value", nobs(fu), "\n")
show("fitted()", fitted(fu))
show("residuals()", residuals(fu))
show("simulate()", simulate(fu, nsim = 1, seed = 1))
show("predict(newdata)", predict(fu, newdata = nd, ndraws = 5))
show("conditional_effects()", conditional_effects(fu))
show("model.frame", model.frame(fu))
show("dharma_residuals", dharma_residuals(fu))
