# Reviewer: the cosy divergences of dev/review-acprior-cosydiv.R with
# the lb = 0 lifted (lb at the window's end, -1/3 for 4 time points,
# which bounds nothing): flat on the natural cosy over the whole window.
# If these have none, the divergences come from the bound at 0.
# Log: dev/review-acprior-log/cosydiv3.txt
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-acprior/dev/review-acprior-cosydiv.R")
cut <- grep("^for [(]s in 3:5[)] [{]", src)[1]
eval(parse(text = src[seq_len(cut - 1L)]))
for (s in 3:5) {
  w <- character(0)
  ds <- withCallingHandlers(
    suppressMessages(frm_sample(fit, chains = 4, iter = 2000, seed = s,
                                refresh = 0, init = "random",
                                prior = set_prior("", class = "cosy",
                                                  lb = -1 / 3))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  dv <- regmatches(w, regexpr("There were [0-9]+ divergent", w))
  cz <- as.matrix(ds)[, "cosy"]
  cat(sprintf("window, seed %d: %s; cosy mean %.4f; share below 0 %.3f\n", s,
              if (length(dv)) dv else "0 divergent", mean(cz), mean(cz < 0)))
}
