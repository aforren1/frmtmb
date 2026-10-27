## brms_fixture(1)'s mo(Exp) frame entry, base against lane: the ported
## brms-suite rows on fixture 1 read it.  Rscript dev/mocode-fixture1.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
d <- as.data.frame(get("brmsfit_example1", envir = asNamespace("brms"))$data)
fr <- suppressWarnings(suppressMessages(frm(bf(count ~ Trt * Age + mo(Exp) +
  s(Age) + volume + offset(Age) + (1 + Trt | visit) + arma(visit, patient),
  sigma ~ Trt), family = student(), data = d, dry_run = "frame")))
mi <- fr$linpreds[[1]]$mo[[1]]
cat(arm, "class", class(d$Exp)[1], "D", mi$D, "codes md5",
    digest::digest(as.integer(mi$codes)), "\n")
