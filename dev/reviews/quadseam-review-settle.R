# Review of lane quadseam: quad_settle()'s acceptance rule read against
# an exact reference. quad_settle() is replaced (in this process only) by
# a copy of the lane's code that also logs each round: the point it left,
# the point the frozen rule's optimizer reached, the re-anchored values it
# compared (f_at, f_new), its decision, and the exact marginal
# log-likelihood (integrate() per group) at both points. A round is a
# wrong decision when it is accepted and the exact objective rose, or
# rejected and the exact objective fell, by more than 1e-9.
# Rounds are forced by a smaller grad_tol (1e-5, 1e-6): the default rarely
# leaves a fit off the re-anchored rule's stationary point.
#   Rscript quadseam-review-settle.R <base|lane> <design,...> <seeds> <tols>
a <- commandArgs(TRUE)
arm <- a[1]
designs <- strsplit(a[2], ",")[[1]]
seeds <- eval(parse(text = a[3]))
tols <- as.numeric(strsplit(a[4], ",")[[1]])
base <- "C:/Users/adf44/source/r/rellib-r20"
lane <- "C:/Users/adf44/source/r/wt-quadseam-lib"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "lane") c(lane, base, ul) else c(base, ul))
suppressMessages(library(frmtmb))
cat("lib:", arm, dirname(find.package("frmtmb")), "\n")
src <- readLines("quadseam-review-acc.R")
# reuse gen(), fo, frm_fam, ref_ll from the accuracy script
eval(parse(text = src[grep("^gen <- function", src):
                        (grep("^for \\(des in designs\\)", src) - 1L)]))
cur <- new.env()
exact_at <- function(par) {
  b <- unname(par[names(par) == "beta"])
  s <- exp(unname(par[names(par) == "theta"][1]))
  -ref_ll(cur$fam, cur$d, b, s)
}
if (arm == "lane") {
  ns <- asNamespace("frmtmb")
  qs <- get("quad_settle", ns)
  body_txt <- deparse(qs)
  # log after each decision: find the line `if (!ok) break`
  i <- grep("^ *if \\(!ok\\) *$", body_txt)
  stopifnot(length(i) == 1L)
  body_txt <- append(body_txt, c(
    "  local({",
    "    e0 <- exact_at(par); e1 <- exact_at(op$par)",
    "    cur$rounds[[length(cur$rounds) + 1L]] <- c(k = k, f_at = f_at,",
    "      f_new = f_new, ok = ok, ex_at = e0, ex_new = e1,",
    "      step = max(abs(op$par - par)))",
    "  })"), after = i - 1L)
  f2 <- eval(parse(text = body_txt))
  environment(f2) <- ns
  # exact_at and cur live in the global env; reach them from the copy
  assign("exact_at", exact_at, envir = globalenv())
  body(f2) <- substitute({
    exact_at <- get("exact_at", envir = globalenv())
    cur <- get("cur", envir = globalenv())
    B
  }, list(B = body(f2)))
  utils::assignInNamespace("quad_settle", f2, ns = "frmtmb")
}
for (des in designs) for (s in seeds) for (tol in tols) {
  g <- gen(des, s)
  cur$fam <- g$fam; cur$d <- g$d; cur$rounds <- list()
  w <- character(); ms <- character()
  t0 <- proc.time()[[3]]
  f <- tryCatch(withCallingHandlers(
    frm(fo[[g$fam]], family = frm_fam[[g$fam]], data = g$d,
        quadrature = TRUE, control = frmtmb_control(grad_tol = tol)),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    },
    message = function(x) {
      ms <<- c(ms, class(x)[1]); invokeRestart("muffleMessage")
    }), error = function(e) e)
  tt <- proc.time()[[3]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("FIT %s %d tol %g %s ERROR %s\n", des, s, tol, arm,
                substr(conditionMessage(f), 1, 200)))
    next
  }
  ex <- exact_at(f$opt$par)
  qsr <- f$opt[["quad_settle"]]
  cat(sprintf(paste0("FIT %s %d tol %g %s obj %.10f exact %.10f d %.3e ",
                     "rounds %s grad %s warn %s msg %s secs %.1f\n"),
              des, s, tol, arm, f$opt$objective, ex, f$opt$objective - ex,
              if (is.null(qsr)) "-" else qsr[["rounds"]],
              if (is.null(qsr)) "-" else format(qsr[["grad"]], digits = 3),
              if (length(w)) paste(w, collapse = ",") else "-",
              if (length(ms)) paste(ms, collapse = ",") else "-", tt))
  for (r in cur$rounds) {
    dex <- r[["ex_new"]] - r[["ex_at"]]
    dq <- r[["f_new"]] - r[["f_at"]]
    wrong <- (r[["ok"]] == 1 && dex > 1e-9) || (r[["ok"]] == 0 && dex < -1e-9)
    cat(sprintf(paste0("  ROUND %d ok %d dq %.3e dexact %.3e step %.3e ",
                       "%s\n"), as.integer(r[["k"]]), as.integer(r[["ok"]]),
                dq, dex, r[["step"]], if (wrong) "WRONG" else ""))
  }
}
cat("END\n")
