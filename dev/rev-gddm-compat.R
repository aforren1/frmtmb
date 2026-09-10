# rev-gddm: is the corrected frm_compat("gddm") predict row accurate,
# or only less wrong? Every clause of the note is run.
#
# Seed 11, the postfit script's own seed and design, so the numbers are
# comparable with dev/gddm-scripts/gddm-postfit-log.txt.
# Arm from GDDM_LIB; default is this review's install of the worktree.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, " frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n\n")

set.seed(11)
ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
n <- 200L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$coh <- factor(rep(c("lo", "hi"), length.out = n))
d$blk <- factor(rep(1:2, each = 2L, length.out = n))
d$cond <- gddm_conditions(d, coh, blk)
fit <- frm(bf(rt | vint(upper, cond) ~ coh, bs ~ blk, ndt ~ 1,
              bias = 0.5),
           family = gddm(control = ctl), data = d)
nd <- unique(d[, c("coh", "blk", "cond")])
nd$rt <- 0.6; nd$upper <- 1L; nd <- nd[order(nd$cond), ]

show <- function(label, e) {
  r <- tryCatch({
    v <- as.numeric(unlist(e))
    paste("ok:", paste(format(head(v, 3), digits = 8), collapse = " "))
  }, error = function(z) paste("FAILS:", conditionMessage(z)))
  cat(sprintf("  %-46s %s\n", label, substr(r, 1, 80)))
}
cat("== every clause of the corrected predict row\n")
show("predict(type = \"link\"), training data",
     predict(fit, type = "link"))
show("predict(type = \"response\"), training data",
     predict(fit, type = "response"))
show("predict(newdata =, type = \"link\")",
     predict(fit, newdata = nd, type = "link"))
show("predict(newdata =, type = \"response\")",
     predict(fit, newdata = nd, type = "response"))
show("conditional_effects()", conditional_effects(fit))
show("fitted()", fitted(fit))
show("residuals(type = \"response\")", residuals(fit))
show("simulate(nsim = 1)", simulate(fit, nsim = 1))
cat("\n  newdata missing a vint() column, the row's other clause:\n")
nd2 <- nd; nd2$cond <- NULL
show("predict(newdata = without cond, type = \"link\")",
     predict(fit, newdata = nd2, type = "link"))
nd3 <- nd; nd3$upper <- NULL
show("predict(newdata = without upper, type = \"link\")",
     predict(fit, newdata = nd3, type = "link"))

cat("\n== the row as frm_compat() prints it\n")
tb <- frm_compat("gddm")
tb <- as.data.frame(tb)
print(tb[tb[[2]] %in% c("predict", "fitted", "simulate", "residuals",
                        "mixture"), 1:3], row.names = FALSE)
