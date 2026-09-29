# Claim 1, the crux: is an UNUSED level of a cs() factor refused with
# the constant-column message, as dev/csfactor-findings.md says, or does
# assemble_frame()'s drop.unused.levels = TRUE remove it first?
#   Rscript dev/csfactor-rev-unused.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
say <- function(...) cat(..., "\n", sep = "")
hdr <- function(x) cat("\n---- ", x, " ----\n", sep = "")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 400L)

set.seed(405)
n <- 240
d <- data.frame(x = rnorm(n))
d$f <- factor(sample(c("a", "b"), n, TRUE), levels = c("a", "b", "c"))
eta <- 0.6 * d$x + c(a = 0, b = 0.8, c = -0.5)[as.character(d$f)]
p1 <- plogis(-0.6 - eta); p2 <- plogis(0.9 - eta)
P <- cbind(p1, p2 - p1, 1 - p2)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))
say("levels(f) = ", paste(levels(d$f), collapse = ","),
    "   table = ", paste(table(d$f), collapse = " "),
    "   level 'c' rows = ", sum(d$f == "c"))

hdr("stats: what model.matrix gives on this factor")
mm <- stats::model.matrix(~ f, d)
say("columns: ", paste(colnames(mm), collapse = ", "))
say("column sums: ", paste(colSums(mm), collapse = ", "))

hdr("frmtmb: does frm() refuse it?")
r <- tryCatch(frm(bf(yo ~ x + cs(f)), family = sratio(), data = d),
              error = function(e) structure(list(m = short(e)),
                                            class = "revfail"))
if (inherits(r, "revfail")) {
  say("REFUSED: ", r$m)
} else {
  say("FITTED. logLik = ", sprintf("%.7f", as.numeric(logLik(r))),
      "  npar = ", length(r$opt$par))
  say("fixef rows: ", paste(rownames(fixef(r)), collapse = ", "))
  say("variables(): ", paste(grep("^bcs|^b_", variables(r), value = TRUE),
                             collapse = ", "))
  dp <- default_prior(r)
  say("default_prior class b coef rows: ",
      paste(dp$coef[dp$class == "b" & nzchar(dp$coef)], collapse = ", "))
  lp <- r$frame$lp[[1L]]
  say("names(frame$lp[[1]]): ", paste(names(lp), collapse = ", "))
  cst <- lp[["cs"]] %||% list()
  say("cs entries in the stored frame = ", length(cst))
  for (ct in cst) {
    say("  entry par=", ct$par, " label=", ct$label,
        " range=", paste(sprintf("%.4g", range(ct$vals)), collapse = ".."),
        " sum=", sprintf("%.4g", sum(ct$vals)))
  }
  mmspec <- lp[["cs_mm"]]
  say("cs_mm terms = ", length(mmspec), "  colnames = ",
      paste(unlist(lapply(mmspec, `[[`, "colnames")), collapse = ", "))
  say("cs_mm xlevels f = ",
      paste(unlist(lapply(mmspec, function(s) s$xlevels$f)),
            collapse = ","))
  hdr("predicting at the unused level 'c'")
  nd <- data.frame(x = 0, f = factor("c", levels = c("a", "b", "c")))
  pr <- tryCatch(fitted(r, newdata = nd),
                 error = function(e) paste0("REFUSED: ", short(e)))
  print(pr)
}

hdr("is the constant-column branch reachable at all?")
d$k <- 5
r2 <- tryCatch(frm(bf(yo ~ x + cs(k)), family = sratio(), data = d),
               error = function(e) paste0("REFUSED: ", short(e)))
say(if (is.character(r2)) r2 else "FITTED (no refusal)")
r3 <- tryCatch(frm(bf(yo ~ cs(k)), family = sratio(), data = d),
               error = function(e) paste0("REFUSED: ", short(e)))
say(if (is.character(r3)) r3 else "FITTED (no refusal)")

hdr("a factor whose unused level is NOT dropped because a level has 1 row")
d2 <- d
d2$f2 <- factor(c(rep("a", n - 1L), "b"), levels = c("a", "b", "c"))
r4 <- tryCatch(frm(bf(yo ~ x + cs(f2)), family = sratio(), data = d2),
               error = function(e) paste0("REFUSED: ", short(e)))
if (is.character(r4)) say(r4) else
  say("FITTED. npar = ", length(r4$opt$par), "  fixef: ",
      paste(rownames(fixef(r4)), collapse = ", "))
cat("\nDONE ", TAG, "\n")
