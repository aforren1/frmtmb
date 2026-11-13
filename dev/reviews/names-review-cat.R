# Review of lane names, 3164: unused categorical levels and K families.
# Usage: Rscript dev/reviews/names-review-cat.R <brms|base|lane>
a <- commandArgs(TRUE)
arm <- a[1]
lib <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-names-lib",
         if (arm != "brms") "C:/Users/adf44/source/r/rellib-r20",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
if (arm == "brms") {
  suppressPackageStartupMessages(library(brms))
} else {
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", find.package("frmtmb"), "\n")
}
set.seed(3164)
n <- 150
d <- data.frame(x = rnorm(n))
eta <- cbind(0, 0.5 + 0.8 * d$x, -0.3 - 0.5 * d$x)
p <- exp(eta) / rowSums(exp(eta))
yy <- apply(p, 1, function(pp) sample(c("a", "b", "c"), 1, prob = pp))
d$y_end <- factor(yy, levels = c("a", "b", "c", "zz"))
d$y_mid <- factor(yy, levels = c("a", "q", "b", "c"))
d$y_first <- factor(yy, levels = c("_0", "a", "b", "c"))
d$y_used <- factor(yy)
d$ynum <- match(yy, c("a", "b", "c")) + 1L       # codes 2..4
d$tr <- 10
d$y1 <- rbinom(n, 10, 0.3)
d$y2 <- rbinom(n, 10 - d$y1, 0.5)
d$y3 <- d$tr - d$y1 - d$y2
d$Yn <- with(d, cbind(y1, y2, y3))
d$Yu <- unname(d$Yn)
q <- function(e) tryCatch(suppressWarnings(suppressMessages(e)),
                          error = function(e) paste("ERROR:",
                                                    conditionMessage(e)))
show <- function(label, v) {
  cat("\n== ", label, "\n", sep = "")
  print(v)
}
if (arm == "brms") {
  dp <- function(f, fam, dd = d) {
    sc <- q(stancode(f, data = dd, family = fam))
    if (startsWith(sc[1], "ERROR")) return(sc)
    sc <- strsplit(sc, "\n")[[1]]
    trimws(grep("vector\\[N\\] mu[^ ]* = ", sc, value = TRUE))
  }
  for (r in c("y_end", "y_mid", "y_first", "y_used", "ynum")) {
    show(paste("categorical", r), dp(as.formula(paste(r, "~ x")),
                                     categorical()))
    sd1 <- q(standata(as.formula(paste(r, "~ x")), data = d,
                      family = categorical()))
    if (is.list(sd1)) show(paste("categorical", r, "Y table"), table(sd1$Y))
  }
  show("multinomial named", dp(bf(Yn | trials(tr) ~ x), multinomial()))
  show("multinomial unnamed", dp(bf(Yu | trials(tr) ~ x), multinomial()))
  fit <- q(brm(y_mid ~ x, data = d, family = categorical(), chains = 1,
               iter = 200, warmup = 100, seed = 1, refresh = 0,
               backend = "rstan"))
  if (!is.character(fit)) {
    show("brms y_mid b_ variables", grep("^b_", variables(fit), value = TRUE))
    pe <- posterior_epred(fit)
    show("brms y_mid epred dimnames", dimnames(pe))
    show("brms y_mid fitted dimnames", dimnames(fitted(fit)))
    pp <- posterior_predict(fit)
    show("brms y_mid posterior_predict values", table(pp))
    show("brms y_mid predict dimnames", dimnames(predict(fit)))
    show("brms y_mid post means", colMeans(as.matrix(fit)[, 1:4]))
  } else show("brms fit", fit)
} else {
  for (r in c("y_end", "y_mid", "y_first", "y_used", "ynum")) {
    f <- q(frm(as.formula(paste(r, "~ x")), data = d,
               family = categorical()))
    if (is.character(f)) { show(paste("categorical", r), f); next }
    show(paste("categorical", r, "fixef"), round(fixef(f)[, 1], 4))
    show(paste("categorical", r, "logLik df"),
         c(ll = as.numeric(logLik(f)), df = attr(logLik(f), "df")))
    show(paste("categorical", r, "fitted dimnames"),
         dimnames(fitted(f))[[3]])
    sim <- q(simulate(f, nsim = 1, seed = 1)[[1]])
    show(paste("categorical", r, "simulate levels"),
         if (is.factor(sim)) levels(sim) else sim)
  }
  for (r in c("y_used", "y_end", "ynum")) {
    f <- q(frm(as.formula(paste(r, "~ x")), data = d,
               family = categorical(K = 3)))
    show(paste("categorical(K = 3)", r),
         if (is.character(f)) f else rownames(fixef(f)))
  }
  f <- q(frm(y_end ~ x, data = d, family = categorical(K = 4)))
  show("categorical(K = 4) y_end",
       if (is.character(f)) f else rownames(fixef(f)))
  for (r in c("Yn", "Yu")) {
    f <- q(frm(as.formula(paste(r, "| trials(tr) ~ x")), data = d,
               family = multinomial()))
    show(paste("multinomial()", r),
         if (is.character(f)) f else rownames(fixef(f)))
    f <- q(frm(as.formula(paste(r, "| trials(tr) ~ x")), data = d,
               family = multinomial(K = 3)))
    show(paste("multinomial(K = 3)", r),
         if (is.character(f)) f else rownames(fixef(f)))
  }
  # a prior written for the K family's dpar names
  f <- q(frm(y_used ~ x, data = d, family = categorical(K = 3),
             prior = set_prior("normal(0, 1)", class = "b", dpar = "mu2")))
  show("categorical(K = 3) with prior on dpar mu2",
       if (is.character(f)) f else rownames(fixef(f)))
  f <- q(frm(bf(y_used ~ x, mu2 ~ 1), data = d, family = categorical(K = 3)))
  show("categorical(K = 3) with mu2 ~ 1 formula",
       if (is.character(f)) f else rownames(fixef(f)))
}
cat("\nDONE\n")
