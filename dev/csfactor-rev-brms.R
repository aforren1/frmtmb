# Claim 4: brms parity without sampling. standata(), get_prior() and
# stancode() only, on the cases the task names: an ordered factor, a
# factor with a level DECLARED but unused, and the newdata-level
# question. Nothing is compiled here.
#   Rscript dev/csfactor-rev-brms.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  # brms FIRST: attaching it after frmtmb masks frmtmb's bf(),
  # default_prior() and the rest
  library(brms)
  library(frmtmb)
})
cat("== build", TAG, ": frmtmb", as.character(packageVersion("frmtmb")),
    " brms", as.character(packageVersion("brms")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 250L)

set.seed(405)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$fc <- factor(sample(c("a", "b", "c"), n, TRUE))
d$fch <- as.character(d$fc)
d$fnum <- factor(as.integer(d$fc))
d$fo <- factor(as.character(d$fc), levels = c("a", "b", "c"),
               ordered = TRUE)
# a level DECLARED but never used
d$fu <- factor(ifelse(d$fc == "c", "b", as.character(d$fc)),
               levels = c("a", "b", "c"))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(d$fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d$yo <- ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L))
cat("table(fu) = ", paste(table(d$fu), collapse = " "),
    "  levels ", paste(levels(d$fu), collapse = ","), "\n", sep = "")

row <- function(lab, bform, ffit) {
  sd <- tryCatch(brms::standata(bform, data = d, family = brms::sratio()),
                 error = function(e) NULL)
  xcs <- if (is.null(sd)) "ERROR" else paste(colnames(sd$Xcs),
                                             collapse = ",")
  xcol <- if (is.null(sd)) "ERROR" else {
    if (is.null(sd$X)) "-" else paste(colnames(sd$X), collapse = ",")
  }
  gp <- tryCatch(brms::get_prior(bform, data = d,
                                 family = brms::sratio()),
                 error = function(e) NULL)
  gpc <- if (is.null(gp)) "ERROR" else {
    paste(gp$coef[gp$class == "b" & nzchar(gp$coef)], collapse = ",")
  }
  f <- tryCatch(suppressWarnings(suppressMessages(ffit())),
                error = function(e) structure(list(m = short(e)),
                                              class = "revfail"))
  ours <- if (inherits(f, "revfail")) paste("REFUSED:", f$m) else {
    paste0(paste(grep("^bcs", variables(f), value = TRUE),
                 collapse = ","),
           "  |fixef| ", nrow(fixef(f)))
  }
  op <- if (inherits(f, "revfail")) "-" else {
    pr <- as.data.frame(default_prior(f))
    paste(pr$coef[pr$class == "b" & nzchar(pr$coef)], collapse = ",")
  }
  cat("\n### ", lab, "\n", sep = "")
  cat("  brms X      : ", xcol, "\n", sep = "")
  cat("  brms Xcs    : ", xcs, "\n", sep = "")
  cat("  brms b coefs: ", gpc, "\n", sep = "")
  cat("  frmtmb      : ", ours, "\n", sep = "")
  cat("  frmtmb b coefs: ", op, "\n", sep = "")
  invisible(f)
}

row("yo ~ x + cs(fc)", brms::bf(yo ~ x + cs(fc)),
    function() frm(frmtmb::bf(yo ~ x + cs(fc)), family = sratio(), data = d))
row("yo ~ x + cs(fch) character", brms::bf(yo ~ x + cs(fch)),
    function() frm(frmtmb::bf(yo ~ x + cs(fch)), family = sratio(), data = d))
row("yo ~ x + cs(fnum)", brms::bf(yo ~ x + cs(fnum)),
    function() frm(frmtmb::bf(yo ~ x + cs(fnum)), family = sratio(), data = d))
row("yo ~ x + cs(fo) ORDERED", brms::bf(yo ~ x + cs(fo)),
    function() frm(frmtmb::bf(yo ~ x + cs(fo)), family = sratio(), data = d))
row("yo ~ x + cs(fu) UNUSED LEVEL", brms::bf(yo ~ x + cs(fu)),
    function() frm(frmtmb::bf(yo ~ x + cs(fu)), family = sratio(), data = d))
row("yo ~ fu + cs(z) unused level in X", brms::bf(yo ~ fu + cs(z)),
    function() frm(frmtmb::bf(yo ~ fu + cs(z)), family = sratio(), data = d))
row("yo ~ cs(fc) + cs(z)", brms::bf(yo ~ cs(fc) + cs(z)),
    function() frm(frmtmb::bf(yo ~ cs(fc) + cs(z)), family = sratio(), data = d))
row("yo ~ x + cs(x) the refused one", brms::bf(yo ~ x + cs(x)),
    function() frm(frmtmb::bf(yo ~ x + cs(x)), family = sratio(), data = d))

cat("\n### ordered factor: are the stored columns brms's Xcs values?\n")
sd <- brms::standata(brms::bf(yo ~ x + cs(fo)), data = d,
                     family = brms::sratio())
fr <- frm(frmtmb::bf(yo ~ x + cs(fo)), family = sratio(), data = d,
          dry_run = "frame")
cs <- fr$linpreds[["yo.mu"]][["cs"]]
cat("  brms Xcs colnames  : ", paste(colnames(sd$Xcs), collapse = ","),
    "\n", sep = "")
cat("  frmtmb cs labels   : ",
    paste(vapply(cs, `[[`, "", "label"), collapse = ","), "\n", sep = "")
for (j in seq_along(cs)) {
  v <- cs[[j]][["vals"]]
  b <- as.numeric(sd$Xcs[, j])
  cat(sprintf("  column %d: max abs difference = %.17g  identical = %s\n",
              j, max(abs(v - b)), identical(v, b)))
}

cat("\n### brms's own Kcs and the stan declaration on the unused level\n")
sdu <- brms::standata(brms::bf(yo ~ x + cs(fu)), data = d,
                      family = brms::sratio())
cat("  brms Kcs = ", sdu$Kcs, "  Xcs colnames = ",
    paste(colnames(sdu$Xcs), collapse = ","), "\n", sep = "")
cat("  brms Xcs column sums = ",
    paste(colSums(sdu$Xcs), collapse = ","), "\n", sep = "")
cat("\nDONE ", TAG, "\n")
