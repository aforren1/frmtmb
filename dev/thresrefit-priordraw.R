## Defect B, both directions, with the RNG pinned so the two arms are
## comparable. FRMTMB_LIB=base selects the 0.64.0 reference build.
arm <- Sys.getenv("FRMTMB_LIB", "lane")
lib <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("ARM", arm, "frmtmb", as.character(packageVersion("frmtmb")), "\n")
say <- function(...) cat(..., "\n", sep = "")
try_msg <- function(expr) {
  tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
}

## direction 1: more than one threshold, so one draw written into the
## whole vector would be the same number repeated
set.seed(701)
db <- data.frame(x = rnorm(40), y = rep(1:4, 10))
for (famnm in c("cratio", "acat", "sratio", "cumulative")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  r <- try_msg(frm_simulate(bf(y ~ x), family = f, data = db,
                            prior = set_prior("normal(0, 2)",
                                              class = "Intercept") +
                              set_prior("normal(0, 1)", class = "b"),
                            nsim = 2, seed = 702))
  say("K=3 ", famnm, " -> ",
      if (is.character(r)) r else "DREW")
  if (!is.character(r)) print(attr(r, "pars"))
}

## direction 2: exactly one threshold. Nothing can recycle, and the
## draw must still happen (the guard must not fail closed)
set.seed(703)
db2 <- data.frame(x = rnorm(40), y = rep(1:2, 20))
for (famnm in c("cratio", "acat", "sratio", "cumulative")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  r <- try_msg(frm_simulate(bf(y ~ x), family = f, data = db2,
                            prior = set_prior("normal(0, 2)",
                                              class = "Intercept") +
                              set_prior("normal(0, 1)", class = "b"),
                            nsim = 2, seed = 704))
  if (is.character(r)) say("K=1 ", famnm, " -> ", r) else {
    p <- attr(r, "pars")
    say("K=1 ", famnm, " -> drew ", paste(names(p), collapse = ","), " = ",
        paste(format(unlist(p[1, ]), digits = 10), collapse = " "),
        "; simulated table = ",
        paste(table(factor(as.integer(r[[1L]]), 1:2)), collapse = "/"))
  }
}

## direction 3: grouped thresholds where EVERY level has one threshold.
## Every entry is then a single index, so the draw runs and each level
## gets its own number
set.seed(705)
dg <- data.frame(x = rnorm(60), g = factor(rep(c("a", "b"), 30)),
                 y = rep(c(1L, 2L, 2L, 1L), 15))
for (famnm in c("cratio", "cumulative")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  r <- try_msg(frm_simulate(bf(y | thres(gr = g) ~ x), family = f,
                            data = dg,
                            prior = set_prior("normal(0, 2)",
                                              class = "Intercept") +
                              set_prior("normal(0, 1)", class = "b"),
                            nsim = 2, seed = 706))
  if (is.character(r)) say("grouped K=1,1 ", famnm, " -> ", r) else {
    p <- attr(r, "pars")
    say("grouped K=1,1 ", famnm, " -> drew ",
        paste(names(p), collapse = ","), " = ",
        paste(format(unlist(p[1, ]), digits = 10), collapse = " "))
  }
}

## direction 4: the internal call, which is where the recycling lives
e <- list(comp = "tau_raw", idx = 1:3, scale = "internal",
          dist = frmtmb:::prior_normal(0, 2))
set.seed(707)
r <- try_msg(frmtmb:::draw_prior_entry(e, "tau_raw[1:3]"))
say("internal draw_prior_entry(idx = 1:3) -> ",
    if (is.character(r)) r else paste("returned one number,",
                                      format(r, digits = 10)))
set.seed(708)
r2 <- try_msg(frmtmb:::draw_prior_pars(list(tau_raw = c(0, 0, 0)),
                                       list(e), "tau_raw[1:3]"))
say("internal draw_prior_pars(idx = 1:3) -> ",
    if (is.character(r2)) r2 else
      paste("tau_raw =", paste(format(r2$est$tau_raw, digits = 10),
                               collapse = " "),
            "| all equal:", length(unique(r2$est$tau_raw)) == 1L))
e1 <- list(comp = "tau_raw", idx = 1L, scale = "internal",
           dist = frmtmb:::prior_normal(0, 2))
set.seed(709)
r3 <- try_msg(frmtmb:::draw_prior_pars(list(tau_raw = 0), list(e1),
                                       "tau_raw_1"))
say("internal draw_prior_pars(idx = 1) -> ",
    if (is.character(r3)) r3 else
      paste("tau_raw =", format(r3$est$tau_raw, digits = 10)))

## and the label, which is what stopped frm_simulate() before the draw
say("prior_entry_label on a 3-index tau_raw entry:")
set.seed(710)
bform <- frmtmb:::as_bform(bf(y ~ x), cratio())
spec <- frmtmb:::parse_spec(bform)
frame <- frmtmb:::assemble_frame(spec, db)
slots <- frmtmb:::nat_slots(frame, list(spec = spec, frame = frame,
                                        bform = bform))
lab <- try_msg(frmtmb:::prior_entry_label(frame, slots, e))
say("  -> ", if (is.character(lab) && length(lab) == 1L) lab else
  paste0("length ", length(lab), ": ", paste(lab, collapse = " ")))
