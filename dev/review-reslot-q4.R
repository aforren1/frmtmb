# Reviewer short check: c() of a set_prior() and a named list. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
dd <- re_data(1L)
m <- tryCatch({ x <- c(set_prior("normal(0, 1)", class = "b"),
                       list(theta_4 = prior_normal(0, 1))); paste("c() gave", class(x)[1]) },
              error = function(e) paste("c() error:", conditionMessage(e)))
cat(m, "\n")
m <- tryCatch({ frm(bf(y ~ x + (1 + x + z | g)), data = dd, prior = x); "fit" },
              error = function(e) paste("frm error:", conditionMessage(e)))
cat(m, "\n")
