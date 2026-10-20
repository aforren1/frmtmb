# Reviewer, re-check of fix round 1: the refit channel cannot leak, the
# new-variable lookup, the MAP refusal, and what update() takes from the
# fit. Arg: lane. Log: dev/drop1raw-rev-out/r1-lane.txt (seed 5, n = 200)
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
ll <- function(f) as.numeric(logLik(f))
tr <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-56s ERROR %s\n", lab,
                substr(gsub("\n", " ", conditionMessage(r)), 1, 170)))
  } else if (inherits(r, "frmtmb_fit")) {
    cat(sprintf("%-56s nobs %d logLik %.10f | call %s\n", lab, nobs(r),
                ll(r), substr(deparse1(getCall(r)), 1, 90)))
  } else {
    cat(sprintf("%-56s %s\n", lab, paste(format(r), collapse = " ")))
  }
  invisible(r)
}
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)))
d$z <- 1 + d$x + log(d$x2) + rnorm(8)[d$g] + rnorm(n)
d$x3[c(3, 17, 40)] <- NA
d$w <- rnorm(n)
d_orig <- d

cat("\n== A. the channel attribute ==\n")
fit <- q(frm(bf(z ~ x + log(x2) + x3), data = d))
tr("fit", fit)
cl <- update(fit, ~ . - x3, evaluate = FALSE)
dd <- cl$data
tr("evaluate = FALSE: class(data), has frmtmb_refit",
   c(class(dd)[1], !is.null(attr(dd, "frmtmb_refit"))))
tr("evaluate = FALSE: nchar(deparse(call))",
   sum(nchar(deparse(cl))))
f_user <- q(frm(bf(z ~ x), data = dd))
tr("user frm(z ~ x, data = <that data frame>)", f_user)
tr("  its na_action rows", paste(unclass(f_user$frame$na_action),
                                 collapse = ","))
tr("  its frame rows", nrow(f_user$frame$data_frame))
u1 <- tr("u1 <- update(fit, ~ . - x3)", update(fit, ~ . - x3))
tr("u1 na_action", paste(unclass(u1$frame$na_action), collapse = ","))
tr("eval(getCall(u1)) in the workspace", eval(getCall(u1)))
tr("u2 <- update(u1, ~ . + x3) (x3 not held by u1)", update(u1, ~ . + x3))
tr("update(u1) (refit of a refit)", update(u1))
tr("update(update(u1, family = student()))",
   update(update(u1, family = student())))
rds <- tempfile(fileext = ".rds"); res <- tempfile(fileext = ".rds")
scr <- tempfile(fileext = ".R")
saveRDS(u1, rds)
writeLines(c(
  paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
  "suppressPackageStartupMessages(library(frmtmb))",
  sprintf("u <- readRDS(%s)", deparse(rds)),
  "out <- list(call = deparse1(getCall(u)), ll = as.numeric(logLik(update(u))),",
  "  ll2 = tryCatch(as.numeric(logLik(update(u, ~ . + x3))),",
  "                 error = function(e) conditionMessage(e)),",
  "  d1 = drop1(u)$AIC)",
  sprintf("saveRDS(out, %s)", deparse(res))), scr)
lg <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(scr)),
              stdout = TRUE, stderr = TRUE)
if (file.exists(res)) {
  o <- readRDS(res)
  cat("fresh session: call", o$call, "\n  update(u) ll", sprintf("%.10f", o$ll),
      "| u1", sprintf("%.10f", ll(u1)), "\n  update(u, ~ . + x3):",
      if (is.character(o$ll2)) substr(o$ll2, 1, 150) else sprintf("%.10f", o$ll2),
      "\n  drop1 AIC", paste(sprintf("%.6f", o$d1), collapse = " "),
      "| here", paste(sprintf("%.6f", q(drop1(u1))$AIC), collapse = " "), "\n")
} else cat("fresh session FAILED:", lg, sep = "\n")

cat("\n== B. new-variable lookup: shadowing objects named d ==\n")
in_fun <- function(fit, dloc) {
  d <- dloc
  tryCatch(q(update(fit, ~ . + w)), error = function(e) e)
}
show_r <- function(lab, r) {
  if (inherits(r, "error")) cat(sprintf("%-56s ERROR %s\n", lab,
    substr(gsub("\n", " ", conditionMessage(r)), 1, 170)))
  else cat(sprintf("%-56s nobs %d logLik %.10f\n", lab, nobs(r), ll(r)))
}
ref <- q(frm(bf(z ~ x + log(x2) + x3 + w), data = d))
show_r("reference frm(... + w) on the global d", ref)
show_r("global d, update(fit, ~ . + w)", tryCatch(q(update(fit, ~ . + w)),
                                                   error = function(e) e))
loc1 <- d; loc1$x2 <- loc1$x2 * 2
show_r("local d: x2 changed", in_fun(fit, loc1))
loc2 <- d; loc2$w <- rnorm(n)
show_r("local d: same model vars, other w", in_fun(fit, loc2))
loc3 <- d[, c("z", "x", "x3", "g", "w")]
show_r("local d: x2 column missing", in_fun(fit, loc3))
x2 <- rexp(n) + 5
show_r("local d: x2 column missing, global x2 vector", in_fun(fit, loc3))
rm(x2)
loc4 <- d[c(2:200, 1), ]
show_r("local d: rows reordered", in_fun(fit, loc4))
loc5 <- d; rownames(loc5) <- paste0("r", 1:n)
show_r("local d: other row names", in_fun(fit, loc5))
loc6 <- d; loc6$x3[c(3, 17, 40)] <- 0; loc6$x3[5] <- NA
show_r("local d: NA moved to another row", in_fun(fit, loc6))
loc7 <- d; loc7$z[3] <- 999
show_r("local d: value changed on a row the fit dropped", in_fun(fit, loc7))
# the fit made inside a function, whose data are a local there: the
# formula environment is that function's frame
mk <- function(dat) {
  dl <- dat
  q(frm(bf(z ~ x + log(x2) + x3), data = dl))
}
fit_l <- mk(d)
show_r("fit with local dl, update(~ . + w) from global",
       tryCatch(q(update(fit_l, ~ . + w)), error = function(e) e))
dl <- d; dl$x2 <- dl$x2 + 1
show_r("  ... with a different global dl", tryCatch(q(update(fit_l, ~ . + w)),
                                                  error = function(e) e))
rm(dl)

cat("\n== C. MAP refusal ==\n")
set.seed(9)
m <- 240
e <- data.frame(x = rnorm(m), x3 = rnorm(m),
                g = factor(sample(sprintf("g%02d", 1:12), m, TRUE)))
e$z <- 1 + e$x + rnorm(12)[e$g] + rnorm(12, 0, 0.5)[e$g] * e$x + rnorm(m)
mp <- function(lab, f, prior) {
  fit <- tryCatch(q(frm(f, data = e, prior = prior)), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-44s FIT ERROR %s\n", lab, conditionMessage(fit)))
    return(invisible())
  }
  f0 <- q(frm(f, data = e))
  r <- tryCatch(q(drop1(fit)), error = function(e) e)
  a <- tryCatch(q(anova(f0, fit)), error = function(e) e)
  cat(sprintf("%-44s logLik-ML %.3e | drop1 %s | anova(ML, it) %s\n", lab,
              ll(fit) - ll(f0),
              if (inherits(r, "error")) "REFUSED" else "runs",
              if (inherits(a, "error")) "REFUSED" else "runs"))
}
mp("sd only, exponential(1)", bf(z ~ x + x3 + (1 | g)),
   set_prior("exponential(1)", class = "sd"))
mp("lkj(2) only", bf(z ~ x + x3 + (1 + x | g)),
   set_prior("lkj(2)", class = "cor"))
mp("lkj(1) only", bf(z ~ x + x3 + (1 + x | g)),
   set_prior("lkj(1)", class = "cor"))
mp("sigma intercept normal(0, 1)", bf(z ~ x + x3, sigma ~ x3),
   set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma"))
mp("bound only, lb = 0 on b_x", bf(z ~ x + x3),
   set_prior("", class = "b", coef = "x", lb = 0))
mp("bound only on sd", bf(z ~ x + x3 + (1 | g)),
   set_prior("", class = "sd", lb = 0))
mp("flat on b written as normal(0, Inf)?", bf(z ~ x + x3),
   tryCatch(set_prior("normal(0, 1e10)", class = "b"), error = function(e) NULL))
# a fit with no prior argument but a mo() term (Dirichlet default?)
e$o <- factor(sample(1:4, m, TRUE), ordered = TRUE)
fmo <- q(frm(bf(z ~ mo(o) + x + x3), data = e))
cat("mo() fit with no prior: prior NULL", is.null(fmo$prior), "| drop1",
    tryCatch({q(drop1(fmo)); "runs"}, error = function(e)
      paste("ERROR", substr(conditionMessage(e), 1, 100))), "\n")
fM <- q(frm(bf(z ~ x + x3), data = e, prior = set_prior("normal(0, 1)",
                                                         class = "b")))
cat("step() on a MAP fit:", tryCatch({q(step(fM, trace = 0)); "runs"},
    error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 90))),
    "\n")
cat("anova(MAP, MAP, refit = TRUE):", tryCatch({q(anova(fM, fM, refit = TRUE));
    "runs"}, error = function(e) "REFUSED"), "\n")

cat("\n== D. what update() takes from the fit ==\n")
fr <- q(frm(bf(z ~ x + x3 + (1 | g)), data = e, REML = TRUE))
tr("update(REML fit)$REML", update(fr)$REML)
tr("update(REML fit, REML = FALSE)$REML", update(fr, REML = FALSE)$REML)
f0 <- q(frm(bf(z ~ x + x3 + (1 | g)), data = e,
            control = frmtmb_control(check_se = "ignore")))
tr("update(fit, control = <other>): check_se of the refit",
   update(f0, control = frmtmb_control(check_se = "stop"))$control$check_se)
tr("update(fit): check_se of the refit", update(f0)$control$check_se)
fb <- q(frm(bf(z ~ x), data = e))
st <- fb$estimates
fs <- local({
  st_l <- st
  q(frm(bf(z ~ x), data = e, start = st_l))
})
tr("start fit: update(fs)", update(fs))
tr("start fit: update(fs, ~ . + x3)", update(fs, ~ . + x3))
tr("reference frm(z ~ x + x3)", frm(bf(z ~ x + x3), data = e))
tr("start fit: update(fs, family = student())",
   update(fs, family = student()))
tr("start fit: update(fs, ~ . + x3, start = NULL)",
   update(fs, ~ . + x3, start = NULL))
