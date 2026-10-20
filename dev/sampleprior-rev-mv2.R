# Reviewer of lane sampleprior: which spelling of a prior on the shared
# nu of a Student-t rescor model the resolver takes.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
f <- bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)
fam <- frmtmb::brmsfamily("student")
for (sp in list(set_prior("gamma(2, 0.1)", class = "nu"),
                set_prior("gamma(2, 0.1)", class = "nu", resp = "y1"),
                set_prior("gamma(2, 0.1)", class = "nu", resp = "y2"))) {
  lab <- paste0("resp='", unclass(sp)[[1]]$resp, "'")
  r <- tryCatch({
    od <- cap_obj(f, D$mv, fam, prior = sp)
    of <- cap_obj(f, D$mv, fam, prior = "flat")
    print(od$env$par)
    p <- od$env$par
    p2 <- p; p2[7] <- p2[7] + 0.5
    sprintf("captured; lp diff %.6f", logprior_at(od, of, p2) -
              logprior_at(od, of, p))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(lab, "->", r, "\n")
}
