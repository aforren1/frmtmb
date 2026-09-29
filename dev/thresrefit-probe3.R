## Probe 3 on the REFERENCE build: settle B, and follow the refit paths
## that reassemble the frame (influence, frm_multiple, update(newdata)).
lib <- Sys.getenv("FRMTMB_PROBE_LIB",
                  "C:/Users/adf44/source/r/rellib-r3")
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("frmtmb", as.character(packageVersion("frmtmb")), "from", lib, "\n")
say <- function(...) cat(..., "\n", sep = "")

## ------------------------------------------------------------------
## B: a class "Intercept" prior draw on an UNORDERED threshold vector.
## The dummy response must show at least two categories, so it cycles.
## ------------------------------------------------------------------
set.seed(1)
ddb <- data.frame(x = rnorm(40), y = rep(1:4, 10))
for (famnm in c("cratio", "acat", "sratio", "cumulative")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  r <- tryCatch({
    s <- frm_simulate(bf(y ~ x), family = f, data = ddb,
                      prior = set_prior("normal(0, 2)",
                                        class = "Intercept") +
                        set_prior("normal(0, 1)", class = "b"),
                      nsim = 3, seed = 5)
    list(ok = TRUE, pars = attr(s, "pars"), y = s)
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  if (!r$ok) {
    say("B: ", famnm, " -> ERROR: ", r$msg)
  } else {
    say("B: ", famnm, " -> pars columns: ",
        paste(names(r$pars), collapse = " | "))
    print(r$pars)
    say("B: ", famnm, " simulated category counts per sim:")
    print(vapply(r$y, function(v) table(factor(as.integer(v),
                                               levels = 1:4)), integer(4)))
  }
}

## The written prior reaches only ONE reported column, so look at what
## the draw actually wrote: re-run the internal pieces directly.
say("---- B: the entry the resolver builds ----")
for (famnm in c("cratio", "cumulative")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  bform <- frmtmb:::as_bform(bf(y ~ x), f)
  spec <- frmtmb:::parse_spec(bform)
  frame <- frmtmb:::assemble_frame(spec, ddb)
  est <- frame[["par_template"]]
  shim <- list(spec = spec, frame = frame, estimates = est)
  pl <- set_prior("normal(0, 2)", class = "Intercept")
  ent <- tryCatch(frmtmb:::resolve_prior_input(shim, pl)$entries,
                  error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(ent)) {
    say("B: ", famnm, " entries -> ", ent)
    next
  }
  for (e in ent) {
    say("B: ", famnm, " entry comp = ", e$comp, " idx = ",
        paste(e$idx, collapse = ","), " scale = ", e$scale,
        " offset = ", if (is.null(e$offset)) "none" else
          paste0(e$offset$comp, "[", paste(e$offset$idx, collapse = ","),
                 "]"))
  }
  set.seed(9)
  dr <- tryCatch(frmtmb:::draw_prior_pars(est, ent, rep("tau", length(ent))),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(dr)) {
    say("B: ", famnm, " draw_prior_pars -> ", dr)
  } else {
    say("B: ", famnm, " tau_raw after the draw = ",
        paste(format(dr$est$tau_raw, digits = 8), collapse = " "))
    say("B: ", famnm, " identical across thresholds = ",
        length(unique(dr$est$tau_raw)) == 1L)
  }
}

## ------------------------------------------------------------------
## A: paths that reassemble the frame
## ------------------------------------------------------------------
mk <- function(seed, n = 40, tau = c(-0.6, 0.5, 2.6), slope = 0.5) {
  set.seed(seed)
  x <- rnorm(n)
  cp <- cbind(plogis(tau[1] - slope * x), plogis(tau[2] - slope * x),
              plogis(tau[3] - slope * x))
  u <- runif(n)
  data.frame(x = x, y = 1L + rowSums(u > cp))
}

## A4: update() with new data that lost the top category
dd <- mk(202)
fit <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = dd))
d2 <- dd
d2$y[d2$y == 4L] <- 3L
say("A4: update(newdata) table = ",
    paste(table(factor(d2$y, levels = 1:4)), collapse = "/"))
u <- tryCatch(suppressWarnings(update(fit, newdata = d2)),
              error = function(e) paste("ERROR:", conditionMessage(e)))
say("A4: update() -> ",
    if (is.character(u)) u else
      paste0("thresholds = ", length(u$estimates$tau_raw)))

## A5: frm_multiple() over imputations, one of which lost the category
imps <- list(dd, d2, dd)
m <- tryCatch(suppressWarnings(
  frm_multiple(bf(y ~ x), data = imps, family = cumulative())),
  error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(m)) {
  say("A5: frm_multiple -> ", m)
} else {
  say("A5: frm_multiple per-fit threshold counts = ",
      paste(vapply(m$fits, function(f) length(f$estimates$tau_raw), 1L),
            collapse = ","))
  s <- tryCatch(summary(m), error = function(e)
    paste("ERROR:", conditionMessage(e)))
  if (is.character(s)) say("A5: summary(frm_multiple) -> ", s) else {
    say("A5: summary rows:")
    print(s)
  }
}

## A6: the autoscale pre-fit
dda <- dd
dda$x <- dda$x * 1e4
fa <- tryCatch(suppressWarnings(
  frm(bf(y ~ x), family = cumulative(), data = dda,
      control = frmtmb_control(autoscale = TRUE))),
  error = function(e) paste("ERROR:", conditionMessage(e)))
say("A6: autoscale fit -> ",
    if (is.character(fa)) fa else
      paste0("thresholds = ", length(fa$estimates$tau_raw)))

## A7: frm_simulate() then frm() -- the user-level two-step
sims <- frm_simulate(bf(y ~ x), family = cumulative(), data = dd,
                     newparams = list(beta = 0.5,
                                      tau_raw = c(-0.6, log(1.1),
                                                  log(2.1))),
                     nsim = 40, seed = 17)
tops <- vapply(sims, function(v) max(as.integer(v)), 1L)
say("A7: frm_simulate seed 17 nsim 40: replicates with max < 4: ",
    sum(tops < 4L))
