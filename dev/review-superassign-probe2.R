# Review of lane superassign: re-entrancy probes. A fit run inside the
# optimizer of another fit (a custom optimizer that fits a second model
# and then runs nlminb on the outer objective), frm_allfit() refits, and
# a nonlinear fit from its default start (nl_fit_off_saddle(),
# nl_restart_round()). Full-precision results to an RDS.
# Usage: Rscript dev/review-superassign-probe2.R <base|lane> <out.rds>
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
R <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "lane") file.path(R, "wt-superassign-lib"),
            file.path(R, "rellib-r18"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
fullp <- function(x) sprintf("%.17g", unlist(x))
res <- list()
set.seed(5)
n <- 300
d <- data.frame(x = runif(n, 0, 3), g = factor(sample(1:15, n, TRUE)))
d$y <- 2 * exp(0.5 * d$x) + rnorm(15, 0, 0.5)[d$g] + rnorm(n, 0, 0.4)
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x))

inner <- list()
nested_opt <- function(par, fn, gr, lower, upper, control) {
  fi <- frm(yc ~ x + (1 | g), family = poisson(), data = d)
  inner[[length(inner) + 1L]] <<- fullp(c(fi$opt$par, fi$opt$objective))
  o <- stats::nlminb(par, fn, gr, lower = lower, upper = upper)
  list(par = o$par, objective = o$objective, convergence = o$convergence)
}
fo <- frm(y ~ x + (1 + x | g), data = d,
          control = frmtmb_control(optimizer = nested_opt))
res$nested <- list(outer = fullp(c(fo$opt$par, fo$opt$objective)),
                   inner = inner)

f0 <- frm(yc ~ x + (1 | g), family = poisson(), data = d)
af <- suppressWarnings(frm_allfit(f0))
res$allfit <- list(print = utils::capture.output(print(af)), num = tryCatch(rapply(lapply(unclass(af), function(z) if (is.list(z) && !is.null(z$opt)) z$opt else z), function(v) sprintf("%.17g", v), classes = "numeric", how = "unlist"), error = function(e) conditionMessage(e)))

# a body on a saddle at its default start: the restart rounds run
msgs <- character(0)
ws <- character(0)
fn <- withCallingHandlers(
  frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), data = d),
  message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  },
  warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
res$nlsaddle <- list(par = fullp(c(fn$opt$par, fn$opt$objective)),
                     msgs = msgs, ws = ws)
saveRDS(res, a[2])
cat("RESULT nested inner fits:", length(inner), "; nl messages:",
    length(msgs), "; nl warnings:", length(ws), "\n")
