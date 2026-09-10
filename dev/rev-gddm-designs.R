# rev-gddm: designs the lane's sweep did not try, and the message a
# user with more than one missing variable is given.
#
# Everything runs at dry_run = "frame". Seed 202609 for the sweep's own
# design, 4242 for the new ones.
# Arm from GDDM_LIB; default is this review's install of the worktree.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, " frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n\n")

ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
run <- function(data, spec, fam = gddm(control = ctl)) {
  tryCatch({
    frm(spec, family = fam, data = data, dry_run = "frame")
    ""
  }, error = function(e) conditionMessage(e))
}
say <- function(label, msg, want) {
  got <- if (nzchar(msg)) "REFUSED" else "accepted"
  cat(sprintf("%s %-50s %s\n", if (got == want) "   " else ">>>",
              label, got))
  if (nzchar(msg))
    cat("      ", substr(gsub("\\s+", " ", msg), 1, 110), "\n", sep = "")
  invisible(msg)
}

set.seed(4242)
n <- 240L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$subj <- factor(sample.int(4L, n, replace = TRUE))
d$coh <- factor(sample(c("lo", "hi"), n, replace = TRUE))
d$blk <- factor(sample.int(2L, n, replace = TRUE))

cat("== A. unbalanced conditions and missing cells\n")
# unbalanced: condition sizes 1, 2, 3, ... by construction
b <- d
b$grp <- factor(rep(seq_len(20L), times = c(1:19, n - sum(1:19))))
b$cond <- gddm_conditions(b, grp)
cat("   condition sizes:", paste(head(as.integer(table(b$cond)), 8),
                                 collapse = " "), "...\n")
say("unbalanced conditions, index correct",
    run(b, bf(rt | vint(upper, cond) ~ grp, bs ~ 1, ndt ~ 1,
              bias = 0.5)), "accepted")
b2 <- b; b2$x <- rnorm(n)
say("unbalanced conditions, x left out of the index",
    run(b2, bf(rt | vint(upper, cond) ~ grp + x, bs ~ 1, ndt ~ 1,
               bias = 0.5)), "REFUSED")

# missing cells: coh x blk with one cell empty
b <- d
keep <- !(b$coh == "hi" & b$blk == "2")
b <- b[keep, ]
b$cond <- gddm_conditions(b, coh, blk)
cat("   cells present:", nrow(unique(b[, c("coh", "blk")])), "of 4,",
    nrow(b), "rows\n")
say("missing cell, index correct",
    run(b, bf(rt | vint(upper, cond) ~ coh * blk, bs ~ 1, ndt ~ 1,
              bias = 0.5)), "accepted")
b3 <- b; b3$cond <- gddm_conditions(b3, coh)
say("missing cell, blk left out of the index",
    run(b3, bf(rt | vint(upper, cond) ~ coh * blk, bs ~ 1, ndt ~ 1,
               bias = 0.5)), "REFUSED")

cat("\n== B. a continuous covariate constant within condition by\n")
cat("      construction, which is the shape a check comparing rows\n")
cat("      rather than values could get wrong\n")
b <- d
b$cond <- gddm_conditions(b, coh, blk)
# one draw per condition, spread back over its rows
u <- stats::runif(length(unique(b$cond)))
b$z <- u[b$cond]
say("z drawn once per condition, mu ~ z",
    run(b, bf(rt | vint(upper, cond) ~ z, bs ~ 1, ndt ~ 1, bias = 0.5)),
    "accepted")
# the same, computed rather than looked up: a group mean is constant
# within its group but is a different double per row
b$zbar <- ave(as.numeric(b$rt), b$cond, FUN = mean)
say("z = ave(rt, cond, mean), a computed constant",
    run(b, bf(rt | vint(upper, cond) ~ zbar, bs ~ 1, ndt ~ 1,
              bias = 0.5)), "accepted")
say("z = ave(rt, cond, mean) but the index is coarser",
    run(within(b, cond <- as.integer(factor(coh))),
        bf(rt | vint(upper, cond) ~ zbar, bs ~ 1, ndt ~ 1, bias = 0.5)),
    "REFUSED")

cat("\n== C. mo(), which is why the check reads columns and not X\n")
b <- d
b$ord <- factor(sample(c("a", "b", "c"), n, replace = TRUE),
                ordered = TRUE)
b$cond <- gddm_conditions(b, ord)
say("mo(ord), index names ord",
    run(b, bf(rt | vint(upper, cond) ~ mo(ord), bs ~ 1, ndt ~ 1,
              bias = 0.5)), "accepted")
b4 <- b; b4$cond <- rep(1:2, length.out = n)
say("mo(ord), ord left out of the index",
    run(b4, bf(rt | vint(upper, cond) ~ mo(ord), bs ~ 1, ndt ~ 1,
               bias = 0.5)), "REFUSED")

cat("\n== D. more than one variable missing from the index\n")
b <- d
b$cond <- rep(1:2, length.out = n)
b$x <- rnorm(n); b$w <- rnorm(n)
m <- say("mu ~ x + w, both missing, bs ~ blk also missing",
         run(b, bf(rt | vint(upper, cond) ~ x + w, bs ~ blk, ndt ~ 1,
                   bias = 0.5)), "REFUSED")
nm <- c("x", "w", "blk")
cat("   variables named in the one message:",
    paste(nm[vapply(nm, function(v)
      grepl(paste0("`", v, "`"), m, fixed = TRUE), TRUE)],
      collapse = " "), "\n")
cat("   parameters named in the one message:",
    paste(c("mu", "bs")[vapply(c("mu", "bs"), function(v)
      grepl(paste0("`", v, "`"), m, fixed = TRUE), TRUE)],
      collapse = " "), "\n")

cat("\n== E. mixture(gddm(), ...), reachability\n")
b <- d; b$cond <- rep(1:2, length.out = n); b$x <- rnorm(n)
for (stage in c("frame", "objective")) {
  r <- tryCatch({
    frm(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1, bias = 0.5),
        family = mixture(gddm(control = ctl), gddm(control = ctl)),
        data = b, dry_run = stage)
    "accepted"
  }, error = function(e) paste("REFUSED:",
                               substr(conditionMessage(e), 1, 90)))
  cat("   dry_run =", stage, ":", r, "\n")
}

cat("\n== F. the other eam families, per row or per group?\n")
# a covariate varying across every grouping must be accepted by each
b <- d; b$cond <- rep(1:2, length.out = n); b$x <- rnorm(n)
b$win <- b$upper + 1L
fams <- list(gaussian = list(gaussian(), quote(rt ~ x + (1 | subj))),
             wiener = list(wiener(),
                           quote(rt | dec(upper) ~ x + (1 | subj))),
             lba = list(lba(2), quote(rt | vint(win) ~ x + (1 | subj))),
             rdm = list(rdm(2), quote(rt | vint(win) ~ x + (1 | subj))),
             wiener_gng = list(wiener_gng(),
                               quote(rt | dec(upper) ~ x + (1 | subj))))
for (nmf in names(fams)) {
  z <- fams[[nmf]]
  r <- tryCatch({
    frm(bf(eval(z[[2]])), family = z[[1]], data = b, dry_run = "frame")
    "accepted"
  }, error = function(e) paste("REFUSED:",
                              substr(conditionMessage(e), 1, 70)))
  cat("   ", nmf, ": ", r, "\n", sep = "")
}
