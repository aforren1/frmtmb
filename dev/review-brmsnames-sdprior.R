## Prior rows beyond class b: sd coefs of transformed group-level terms,
## a nonlinear parameter's transformed term, category-specific cs().
## brms 2.23.0 default_prior() (no compile) against frmtmb's, and every
## brms coef row (and the old design-column spelling) through
## frmtmb's validate_prior(). Seed 9.
##   Rscript dev/review-brmsnames-sdprior.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(9)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:12, 20)))
d$y <- 1 + d$x + 0.3 * d$z^2 + rnorm(12)[d$g] + rnorm(n)
d$yo <- cut(d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
cases <- list(
  list(f = "y ~ I(x^2) + (1 + I(x^2) | g)", fam = "gaussian",
       old = list(class = "sd", coef = "I(x^2)", group = "g")),
  list(f = "y ~ poly(x, 2) + (1 + poly(x, 2) | g)", fam = "gaussian",
       old = list(class = "sd", coef = "poly(x, 2)2", group = "g")),
  list(f = "y ~ x + (1 + exp(x) || g)", fam = "gaussian",
       old = list(class = "sd", coef = "exp(x)", group = "g")),
  list(f = "bf(y ~ a + b * x, a ~ I(z^2), b ~ 1, nl = TRUE)",
       fam = "gaussian",
       old = list(class = "b", coef = "I(z^2)", nlpar = "a")),
  list(f = "bf(yo ~ cs(I(x^2)) + z)", fam = "acat",
       old = list(class = "b", coef = "I(x^2)"))
)
key <- function(p) with(p, paste(class, coef, group, resp, dpar, nlpar,
                                 sep = "|"))
brms_tabs <- lapply(cases, function(cs) {
  f <- eval(parse(text = sub("^bf\\(", "brms::bf(", cs$f)))
  if (!inherits(f, "brmsformula")) f <- brms::bf(f)
  as.data.frame(q(brms::default_prior(f, data = d, family = cs$fam)))
})
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
for (i in seq_along(cases)) {
  cs <- cases[[i]]
  f <- eval(parse(text = cs$f))
  if (!inherits(f, "brmsformula")) f <- bf(f)
  fam <- get(cs$fam, asNamespace("frmtmb"))()
  bp <- brms_tabs[[i]]
  fp <- tryCatch(as.data.frame(q(default_prior(f, data = d, family = fam))),
                 error = function(e) e)
  cat("==", cs$f, "\n")
  if (inherits(fp, "error")) {
    cat("  frmtmb default_prior ERROR:", conditionMessage(fp), "\n")
    next
  }
  kb <- key(bp)
  kf <- key(fp)
  cat("  rows only brms:", setdiff(kb, kf), "\n")
  cat("  rows only frm: ", setdiff(kf, kb), "\n")
  for (r in which(nzchar(bp$coef))) {
    a <- list(prior = if (bp$class[r] == "sd") "student_t(3, 0, 5)" else
      "normal(0, 5)", class = bp$class[r], coef = bp$coef[r])
    for (nm in c("group", "dpar", "nlpar", "resp")) {
      if (nzchar(bp[[nm]][r])) a[[nm]] <- bp[[nm]][r]
    }
    v <- tryCatch(q(validate_prior(do.call(set_prior, a), f, data = d,
                                   family = fam)), error = function(e) e)
    res <- if (inherits(v, "error")) conditionMessage(v) else {
      v <- as.data.frame(v)
      if (any(v$source == "user" & v$coef == a$coef)) "ok" else
        "not applied"
    }
    if (res != "ok") cat("  brms coef", a$class, a$coef, "->", res, "\n")
  }
  a <- c(list(prior = "normal(0, 5)"), cs$old)
  if (a$class == "sd") a$prior <- "student_t(3, 0, 5)"
  v <- tryCatch(q(validate_prior(do.call(set_prior, a), f, data = d,
                                 family = fam)), error = function(e) e)
  cat("  old spelling", a$class, shQuote(a$coef), "->",
      if (inherits(v, "error")) conditionMessage(v) else {
        v <- as.data.frame(v)
        u <- v[v$source == "user", ]
        paste(nrow(v), "rows; user rows", paste(u$coef, collapse = ","))
      }, "\n")
}
