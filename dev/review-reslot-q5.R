# Reviewer short check of fix round 2: the route mark. Usage:
# Rscript dev/review-reslot-q5.R brms|after. Seed 1.
arm <- commandArgs(trailingOnly = TRUE)[1L]
set.seed(1)
G <- 10L; n <- 8L
d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[d$g] + rnorm(G * n)
f1 <- y ~ x + (1 + x | g)
out <- "C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/review-reslot-out/r3"
dir.create(out, showWarnings = FALSE)
if (identical(arm, "brms")) {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r14",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages(library(brms))
  tb <- default_prior(bf(f1), data = d)
  tb$prior[tb$class == "b" & tb$coef == "x"] <- "normal(0, 1)"
  saveRDS(tb, file.path(out, "brms-table.rds"))
  cat("saved brms table\n")
  quit(save = "no")
}
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
q <- function(expr) suppressWarnings(suppressMessages(expr))
res <- function(expr) {
  r <- tryCatch(q(expr), error = function(e) e)
  if (inherits(r, "error")) {
    m <- conditionMessage(r)
    if (grepl("route = \"fit\"", m, fixed = TRUE)) "REFUSED (fit-route)" else
      paste("ERROR:", substr(m, 1, 160))
  } else r
}
ents <- function(rp) {
  if (!is.list(rp) || is.null(rp$ri)) return(if (is.character(rp)) rp else "none")
  paste(vapply(rp$ri$entries, function(e) {
    paste0(e$comp, "[", paste(e$idx, collapse = ","), "]:", e$dist$kind,
           if (isTRUE(e$dist$flat)) "-flat" else "")
  }, ""), collapse = " ")
}
fit <- q(frm(f1, data = d))
rs <- function(pr) ents(res(frmtmb.sample:::sample_resolve_priors(fit, pr)))
p <- default_prior(f1, data = d)
pe <- p; pe$prior[pe$class == "b" & pe$coef == "x"] <- "normal(0, 1)"
cat("== S1: frm()'s table on frm_sample, every edit\n")
cases <- list(
  "unedited" = p,
  "$<- one row" = pe,
  "[<- one cell" = { z <- p; z[2, "prior"] <- "normal(0, 1)"; z },
  "subset to the edited row" = pe[pe$class == "b" & pe$coef == "x", ],
  "subset to rows with a density" = pe[pe$prior != "(flat)" & pe$prior != "", ],
  "rbind(pe, pe[1, ])" = rbind(pe, pe[1, ]),
  "as.brmsprior(pe)" = as.brmsprior(pe),
  "set_prior() + as.brmsprior(pe)" = set_prior("normal(0, 2)", class = "sd") +
    as.brmsprior(pe),
  "c(as.brmsprior(pe), set_prior())" = res(c(as.brmsprior(pe),
                                              set_prior("normal(0, 2)",
                                                        class = "sd"))),
  "validate_prior(sp, f) (fit)" = validate_prior(set_prior("normal(0, 1)",
                                                           class = "b"),
                                                 f1, data = d),
  "transform(pe) (plain data.frame?)" = transform(pe),
  "pe[order(pe$class), ]" = pe[order(pe$class), ],
  "within/merge-free: pe with a column added" = { z <- pe; z$note <- "a"; z },
  "as.brmsprior(as.data.frame(pe)) (unmarked)" = as.brmsprior(as.data.frame(pe)))
for (nm in names(cases)) {
  x <- cases[[nm]]
  cat(sprintf("  %-42s route=%-6s -> %s\n", nm,
              format(attr(x, "route") %||% "<none>"),
              if (is.character(x)) x else rs(x)))
}
cat("\n== the frm_sample() entry points\n")
cat("  frm_sample(fit, prior = pe): ",
    substr(format(res(frm_sample(fit, prior = pe, chains = 1, iter = 20,
                                 refresh = 0))), 1, 60)[1], "\n")
cat("  frm_sample(formula, prior = pe): ",
    substr(format(res(frm_sample(f1, data = d, prior = pe, chains = 1,
                                 iter = 20, refresh = 0))), 1, 60)[1], "\n")
ds <- q(frm_sample(fit, chains = 1, iter = 60, warmup = 30, seed = 1,
                   refresh = 0))
cat("  update(draws, prior = pe): ",
    substr(format(res(update(ds, prior = pe, chains = 1, iter = 20,
                             refresh = 0))), 1, 60)[1], "\n")
cat("  route of the run's prior_summary():",
    format(attr(prior_summary(ds), "route") %||% "<none>"), "\n")
cat("  prior_summary() passed back:", rs(prior_summary(ds)), "\n")
cat("  defaults for comparison:    ", rs(NULL), "\n")

cat("\n== the mark and frm()\n")
same <- function(a, b) identical(a$opt$par, b$opt$par) &&
  identical(a$opt$objective, b$opt$objective)
f_pe <- q(frm(f1, data = d, prior = pe))
f_sp <- q(frm(f1, data = d, prior = set_prior("normal(0, 1)", class = "b",
                                               coef = "x")))
cat("  frm(prior = pe) == frm(prior = b x spec):", same(f_pe, f_sp), "\n")
cat("  route on fit$prior:", format(attr(f_pe$prior, "route") %||% "<none>"),
    "\n")
cat("  frm_sample(MAP fit made from pe), own prior: ",
    rs(NULL) == ents(res(frmtmb.sample:::sample_resolve_priors(
      f_pe, NULL, base = f_pe$prior))), "\n")
cat("    entries:", ents(res(frmtmb.sample:::sample_resolve_priors(
  f_pe, NULL, base = f_pe$prior))), "\n")
f_ml <- q(frm(f1, data = d))
f_p <- q(frm(f1, data = d, prior = p))
cat("  frm(prior = unedited fit table) == frm():", same(f_ml, f_p), "\n")
cat("  update(frm fit, prior = pe) == frm(prior = b x spec):",
    same(q(update(f_ml, prior = pe)), f_sp), "\n")
ps <- default_prior(f1, data = d, route = "sample")
f_s <- q(frm(f1, data = d, prior = ps))
cat("  frm(prior = sample-route table): fits, prior classes:",
    paste(unique(as.data.frame(f_s$prior)$class), collapse = ","), "\n")

cat("\n== tables read as written\n")
pse <- ps; pse$prior[pse$class == "b" & pse$coef == "x"] <- "normal(0, 1)"
cat("  sample-route table, edited:", rs(pse), "\n")
bt <- readRDS(file.path(out, "brms-table.rds"))
cat("  brms default_prior(), edited:", rs(bt), "\n")
hb <- as.brmsprior(data.frame(prior = c("normal(0, 1)", ""), class = c("b", "sd"),
                              coef = c("x", ""), stringsAsFactors = FALSE))
cat("  hand-built (b x normal, sd ''):", rs(hb), "\n")

cat("\n== S2: c() and list() mixed inputs\n")
m <- function(expr) tryCatch({ force(expr); "ACCEPTED" },
                             error = function(e) conditionMessage(e))
cat("  c(sp, list(theta_4 = )):",
    m(c(set_prior("normal(0, 1)", class = "b"),
        list(theta_4 = prior_normal(0, 1)))), "\n")
cat("  c(sp, prior_normal()):",
    m(c(set_prior("normal(0, 1)", class = "b"), prior_normal(0, 1))), "\n")
cat("  frm(prior = list(sp, theta_4 = )):",
    m(q(frm(f1, data = d, prior = list(set_prior("normal(0, 1)", class = "b"),
                                       theta_4 = prior_normal(0, 1))))), "\n")
cat("  frm(prior = list(sp, sp2)) (unnamed list of specs):",
    m(q(frm(f1, data = d, prior = list(set_prior("normal(0, 1)", class = "b"),
                                       set_prior("normal(0, 2)",
                                                 class = "sd"))))), "\n")
