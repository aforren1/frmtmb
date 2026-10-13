# Reviewer, re-check: frmtmb's update(prior = ) on the fixed build, the
# same four attacks as dev/moparity-rev2-brms-update.R.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(11)
dd <- data.frame(x = rnorm(60), z = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
p0 <- set_prior("normal(0, 0.05)", class = "b", coef = "x") +
  set_prior("normal(0, 10)", class = "Intercept")
f0 <- frm(bf(y ~ x + z), family = gaussian(), data = dd, prior = p0)
rows <- function(f) {
  if (is.null(f$prior)) return("<none>")
  d <- as.data.frame(f$prior)
  paste(paste0(d$class, "/", d$coef, "/", d$dpar, ":", d$prior), collapse = "; ")
}
bx <- function(f) if ("x" %in% rownames(fixef(f))) format(fixef(f)["x", 1], digits = 6) else "(no x)"
try1 <- function(lab, expr) {
  msgs <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  if (inherits(r, "error")) cat("==", lab, ": ERROR", substr(conditionMessage(r), 1, 160), "\n")
  else cat("==", lab, ": b_x", bx(r), "| prior", rows(r), "\n")
  if (length(msgs)) cat("   message:", substr(paste(msgs, collapse = " | "), 1, 200), "\n")
}
cat("stored: b_x", bx(f0), "| prior", rows(f0), "\n")
# frmtmb's own table, as get_prior() gives it
tbl <- get_prior(bf(y ~ x + z), family = gaussian(), data = dd)
cat("class of get_prior():", class(tbl)[1], "| priors in it:",
    paste(unique(tbl$prior), collapse = ", "), "\n")
t1 <- tbl
t1$prior[t1$class == "sigma"] <- "exponential(2)"
try1("frmtmb table, sigma edited", update(f0, prior = t1))
t2 <- tbl[tbl$class == "b" & tbl$coef == "x", ]
t2$prior <- ""
try1("frmtmb table, one row b x, prior empty", update(f0, prior = t2))
try1("frmtmb whole table unedited", update(f0, prior = tbl))
# brms's own table
bt <- brms::default_prior(y ~ x + z, data = dd)
cat("brms table priors:", paste(paste0(bt$class, "/", bt$coef, ":", bt$prior), collapse = "; "), "\n")
bt1 <- bt
bt1$prior[bt1$class == "sigma"] <- "exponential(2)"
try1("brms table, sigma edited", update(f0, prior = bt1))
bt2 <- bt[bt$class == "b" & bt$coef == "x", ]
try1("brms table, one row b x, prior empty", update(f0, prior = bt2))
try1("duplicate slots in new", update(f0, prior =
  set_prior("normal(0, 1)", class = "b", coef = "x") +
  set_prior("normal(0, 2)", class = "b", coef = "x")))
try1("duplicate slots in new (brms)", update(f0, prior =
  brms::prior(normal(0, 1), class = b, coef = x) +
  brms::prior(normal(0, 2), class = b, coef = x)))
try1("formula drops x, new row on b x", update(f0, formula. = ~ . - x,
  prior = set_prior("normal(0, 1)", class = "b", coef = "x")))
try1("formula drops x, new sigma row", update(f0, formula. = ~ . - x,
  prior = set_prior("exponential(2)", class = "sigma")))
try1("empty brmsprior", update(f0, prior = brms::empty_prior()))
