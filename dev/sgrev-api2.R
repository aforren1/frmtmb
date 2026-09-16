# Leftovers: the rendered Rd section read in full, bbmle in the order
# the lane's defect 4 names, and whether any of the ten positional
# divergences answers a different question SILENTLY instead of
# erroring.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))

cat("== the rendered borrowed-generic seam, in full ==\n")
q(library(frmtmb))
rd <- tools::Rd_db("frmtmb")[["frmtmb-sampling-api.Rd"]]
f <- tempfile()
tools::Rd2txt(rd, out = f, options = list(underline_titles = FALSE))
ln <- readLines(f, warn = FALSE)
i <- grep("borrowed-generic seam", ln)
j <- grep("^[A-Za-z].*:$", ln)
j <- j[j > i[1]]
end <- if (length(j)) j[1] - 1L else min(i[1] + 40L, length(ln))
cat(paste(ln[i[1]:end], collapse = "\n"), "\n")

cat("\n== bbmle, in the order defect 4 names ==\n")
q(library(bbmle)); q(library(frmtmb.sample))
cat("search head: ", paste(head(search(), 4), collapse = " "), "\n")
g <- get("parnames", envir = globalenv())
cat("parnames resolves to: ", environmentName(environment(g)),
    "   formals ", paste(deparse(args(g)), collapse = " "), "\n")
x <- numeric(3)
cat("plain parnames(x) (bbmle's meaning is 'give me the names'): ",
    tryCatch(paste(capture.output(print(parnames(x))), collapse = ""),
             error = function(e) paste("ERROR:",
               substr(conditionMessage(e), 1, 60))), "\n")
cat("bbmle::parnames(x): ",
    tryCatch(paste(capture.output(print(bbmle::parnames(x))),
                   collapse = ""),
             error = function(e) paste("ERROR:",
               substr(conditionMessage(e), 1, 60))), "\n")
parnames(x) <- c("a", "b", "c")
cat("parnames(x) <- ... still works (the replacement name is not",
    "masked): ", !is.null(attr(x, "parnames")), "\n")
cat("and then plain parnames(x) reads: ",
    tryCatch(paste(capture.output(print(parnames(x))), collapse = ""),
             error = function(e) paste("ERROR:",
               substr(conditionMessage(e), 1, 60))), "\n")
cat("bbmle::parnames(x) reads: ",
    tryCatch(paste(capture.output(print(bbmle::parnames(x))),
                   collapse = ""),
             error = function(e) "ERROR"), "\n")

cat("\n== the ten positional divergences: error or silence? ==\n")
q(library(frmtmb))
ds <- readRDS("dev/stan-cache/sgrev-draws.rds")
nd <- ds$fit$frame$data[1:5, , drop = FALSE]
try1 <- function(lbl, e) {
  r <- tryCatch(suppressWarnings(suppressMessages(eval(e))),
                error = function(x) structure(conditionMessage(x),
                                              class = "e"))
  cat(sprintf("  %-40s %s\n", lbl,
      if (inherits(r, "e")) paste("ERROR:", substr(unclass(r), 1, 46))
      else paste("RETURNED", class(r)[1], "of length", length(r))))
}
try1("as.mcmc(ds, TRUE)", quote(as.mcmc(ds, TRUE)))
try1("log_lik(ds, nd)", quote(log_lik(ds, nd)))
try1("mcmc_plot(ds, 'b_x')", quote(mcmc_plot(ds, "b_x")))
try1("posterior_epred(ds, nd, NA)", quote(posterior_epred(ds, nd, NA)))
try1("posterior_interval(ds, 0.9)", quote(posterior_interval(ds, 0.9)))
try1("posterior_linpred(ds, TRUE, nd, NA)",
     quote(posterior_linpred(ds, TRUE, nd, NA)))
try1("posterior_predict(ds, nd, NA)", quote(posterior_predict(ds, nd, NA)))
try1("pp_mixture(ds, nd)", quote(pp_mixture(ds, nd)))
try1("predictive_error(ds, nd)", quote(predictive_error(ds, nd)))
try1("psis(ds, nd)", quote(psis(ds, nd)))
cat("DONE\n")
