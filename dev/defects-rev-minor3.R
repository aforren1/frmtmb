# Reviewer of lane defects, recheck of minor 3: default_prior() and
# validate_prior() on responses the family cannot read, frmtmb (lane)
# against brms 2.23.0.
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({loadNamespace("frmtmb"); loadNamespace("brms")})
set.seed(111)
d <- data.frame(x = rnorm(30), yb = rnorm(30), y3 = sample(0:2, 30, TRUE),
                cnt = rnorm(30), yc = rnorm(30), yord = sample(1:3, 30, TRUE))
cases <- list(
  beta = list(yb ~ x, "Beta"), bern = list(y3 ~ x, "bernoulli"),
  pois = list(cnt ~ x, "poisson"), cumul_cont = list(yc ~ x, "cumulative"),
  sratio_cont = list(yc ~ x, "sratio"), cumul_ok = list(yord ~ x, "cumulative"))
key <- function(t) {
  t <- as.data.frame(t)
  apply(t[c("class", "coef", "group", "resp", "dpar", "nlpar")], 1,
        function(r) paste(ifelse(is.na(r), "", r), collapse = "|"))
}
m <- function(e) tryCatch(e, error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 110)))
for (nm in names(cases)) {
  fo <- cases[[nm]][[1]]; fam <- cases[[nm]][[2]]
  bfam <- get(fam, asNamespace("brms"))()
  ffam <- if (fam == "poisson") stats::poisson() else get(fam, asNamespace("frmtmb"))()
  bd <- m(key(suppressWarnings(brms::default_prior(fo, data = d, family = bfam))))
  fd <- m(key(suppressWarnings(frmtmb::default_prior(fo, data = d, family = ffam))))
  bv <- m(key(suppressWarnings(brms::validate_prior(brms::set_prior("normal(0,1)", class = "b"), fo, data = d, family = bfam))))
  fv <- m(key(suppressWarnings(frmtmb::validate_prior(frmtmb::set_prior("normal(0,1)", class = "b"), fo, data = d, family = ffam))))
  st <- function(v) if (length(v) == 1 && grepl("^ERROR", v)) v else paste(length(v), "rows")
  cat(sprintf("%-12s default_prior brms: %s | frmtmb: %s\n", nm, st(bd), st(fd)))
  cat(sprintf("%-12s validate_prior brms: %s | frmtmb: %s\n", "", st(bv), st(fv)))
  if (!grepl("^ERROR", bd[1]) && !grepl("^ERROR", fd[1])) {
    cat("             common rows:", sum(fd %in% bd), " only brms:", paste(setdiff(bd, fd), collapse = " ; "), "\n")
  }
  cat("             frm():", substr(m(class(suppressWarnings(frmtmb::frm(fo, data = d, family = ffam)))), 1, 90), "\n")
}
