# Lane optrap, item 4: after an optimizer run returns, is the objective's
# last point (where parList() reads the random effects) and best point
# (where sdreport() reads) the optimum it reports? Traced at the exit of
# optimize_obj() and of fit_assembled() on a few model families.
#   Rscript dev/optrap-state.R <lib or "">
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
at <- function(obj, par) {
  lf <- obj$env$lfixed()
  c(last = identical(unname(obj$env$last.par[lf]), unname(par)),
    best = identical(unname(obj$env$last.par.best[lf]), unname(par)))
}
tally <- new.env()
tally$opt <- list()
trace("optimize_obj", where = ns, print = FALSE,
      exit = quote({
        r <- returnValue()
        if (!is.null(r$par) && length(r$par)) {
          k <- length(get("opt", envir = .GlobalEnv$tally)) + 1L
          assign("opt", c(get("opt", envir = .GlobalEnv$tally),
                          list(.GlobalEnv$at(obj, r$par))),
                 envir = .GlobalEnv$tally)
        }
      }))
fits <- list()
data(sleepstudy, package = "lme4")
cases <- list(
  sleep = function() frm(Reaction ~ Days + (Days | Subject),
                         data = sleepstudy),
  pois = function() {
    set.seed(2)
    d <- data.frame(g = factor(rep(1:30, each = 8)), x = rnorm(240))
    d$y <- rpois(240, exp(0.5 + 0.3 * d$x + rnorm(30, 0, 0.6)[d$g]))
    frm(y ~ x + (1 | g), family = poisson(), data = d)
  },
  binom = function() {
    set.seed(3)
    d <- data.frame(g = factor(rep(1:30, each = 8)), x = rnorm(240))
    d$y <- rbinom(240, 1, plogis(-0.2 + 0.8 * d$x + rnorm(30, 0, 1)[d$g]))
    frm(y ~ x + (1 | g), family = bernoulli(), data = d)
  },
  ri20 = function() {
    set.seed(19)
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + 0.5 * d$x + rnorm(100)
    frm(y ~ x + (1 | g), data = d)
  },
  skewn = function() {
    set.seed(20)
    g <- factor(rep(1:15, each = 10))
    d <- data.frame(g = g, x = rnorm(150))
    d$y <- 1 + 0.5 * d$x + rnorm(15, sd = 0.8)[g] + rnorm(150)
    frm(y ~ x + (1 | g), family = skew_normal(), data = d)
  })
for (nm in names(cases)) {
  tally$opt <- list()
  f <- suppressWarnings(suppressMessages(cases[[nm]]()))
  st <- do.call(rbind, tally$opt)
  fin <- at(f$obj, f$opt$par)
  ref <- ns$solved_par_list(f$obj, f$opt$par)
  cat(sprintf(paste0("%-6s runs %d | at exit of optimize_obj: last %d of %d,",
                     " best %d of %d | fit: last %s best %s | ",
                     "max |ranef - modes at the optimum| %.3g\n"),
              nm, nrow(st), sum(st[, "last"]), nrow(st), sum(st[, "best"]),
              nrow(st), fin[["last"]], fin[["best"]],
              max(abs(unname(f$estimates$b) - unname(ref$b)))))
}
