# Reviewer of lane sampleprior: for each reviewer case, the parameter
# vector frm_sample() samples and brms's Stan lprior block.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
which <- commandArgs(TRUE)
if (!length(which)) which <- names(CASES)
for (nm in which) {
  cs <- CASES[[nm]]
  f <- eval(cs$frm)
  b <- eval(cs$brm)
  cat("\n=====", nm, "=====\n")
  od <- tryCatch(cap_obj(f[[1]], f[[2]], f[[3]]),
                 error = function(e) conditionMessage(e))
  if (is.character(od)) { cat("FRM ERROR:", od, "\n") } else {
    print(od$env$par)
    cat("random:", od$env$random, "\n")
  }
  sc <- tryCatch(brms::stancode(b[[1]], data = f[[2]], family = b[[2]]),
                 error = function(e) conditionMessage(e))
  s <- strsplit(as.character(sc), "\n")[[1]]
  i1 <- grep("^parameters \\{", s)
  i2 <- grep("^model \\{", s)
  if (length(i1) && length(i2)) cat(s[i1:(i2 - 1)], sep = "\n")
  else cat(s, sep = "\n")
}
