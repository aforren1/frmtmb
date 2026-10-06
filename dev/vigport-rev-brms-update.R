# Reviewer: defect 4. What does brms 2.23.0's update() do with a
# class-wide lkj(2) prior when the new formula drops the last
# correlation? A real compile and a short run (1 chain of 300), since
# update() of an empty = TRUE fit errors in brms; the lane's brms-vig
# run (overview.7.1 OK, Rhat 1.003) already shows the refit samples.
#
#   Rscript dev/vigport-rev-brms-update.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
data("kidney", package = "brms")
pr <- c(set_prior("normal(0,5)", class = "b"),
        set_prior("cauchy(0,2)", class = "sd"),
        set_prior("lkj(2)", class = "cor"))
f1 <- brm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal(), prior = pr, chains = 1,
          iter = 300, refresh = 0, seed = 1)
cat("\nfit1 user priors:\n")
print(f1$prior[f1$prior$source == "user", c("prior", "class", "coef",
                                            "group")])
f2 <- tryCatch(
  withCallingHandlers(
    update(f1, formula. = ~ . - (1 + age | patient) + (1 | patient),
           chains = 1, iter = 300, refresh = 0, seed = 1),
    warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }, message = function(m) {
      cat("MESSAGE:", conditionMessage(m))
      invokeRestart("muffleMessage")
    }),
  error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  })
if (!is.null(f2)) {
  cat("\nupdate() returned a", class(f2)[1], "\nfit2 user priors:\n")
  print(f2$prior[f2$prior$source == "user", c("prior", "class", "coef",
                                              "group")])
  cat("fit2 formula:", deparse1(formula(f2)$formula), "\n")
  cat("lkj in fit2's Stan code:", grepl("lkj", stancode(f2)), "\n")
}
cat("\n## frmtmb, the same update\n")
suppressMessages(library(frmtmb))
g1 <- frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal(), prior = pr)
r <- tryCatch({
  update(g1, formula. = ~ . - (1 + age | patient) + (1 | patient))
  "ACCEPTED"
}, error = function(e) paste("REFUSED:", conditionMessage(e)))
cat(r, "\n")
# the guard's absent case: drop only the class-wide lkj and update
pr2 <- pr[pr$class != "cor", ]
g2 <- frm(time | cens(censored) ~ age * sex + disease + (1 + age | patient),
          data = kidney, family = lognormal(), prior = pr2)
r2 <- tryCatch({
  update(g2, formula. = ~ . - (1 + age | patient) + (1 | patient))
  "ACCEPTED"
}, error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("without the lkj prior:", r2, "\n")
