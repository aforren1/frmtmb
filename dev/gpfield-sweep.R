# Lane gpfield: the lane against 0.82.0's route, interleaved row by row in
# one process on wt-gpfield-lib, the order alternating by row. Arm "base"
# mocks this lane's changes off: gp_mean_relin_applies() and
# gp_hsgp_relin() return nothing (no second reading of a mean or of an
# approximate field) and gp_end_messages() returns no warning. The kernel
# code is not on these rows' path (exp_quad keeps 0.82.0's arithmetic).
# Args: set part
#   set sweep122: lane gpstop's review sweep (dev/gpstop-fx2-sweep.R mk()),
#     part gauss (gaussian n 40, 120, seeds 1, 2, and the 7 k = 12 rows), pb
#     (poisson, bernoulli), student (student, and the 4 n = 400 rows)
#   set fresh, lane96: fitedge's 192 rows (dev/gpstop-gpfresh.R), part all
#   set ng150: lane arstep's 150 rows (dev/gpstop-ng.R, seed 1), part a
#     comma list of cases
a <- commandArgs(TRUE)
set <- a[1]; part <- a[2]
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
ns <- asNamespace("frmtmb")
real <- list(applies = ns$gp_mean_relin_applies, hsgp = ns$gp_hsgp_relin,
             msgs = ns$gp_end_messages)
arm_set <- function(arm) {
  if (arm == "base") {
    gf_put("gp_mean_relin_applies", function(...) FALSE)
    gf_put("gp_hsgp_relin", function(...) NULL)
    gf_put("gp_end_messages", function(fit) list())
  } else {
    gf_put("gp_mean_relin_applies", real$applies)
    gf_put("gp_hsgp_relin", real$hsgp)
    gf_put("gp_end_messages", real$msgs)
  }
}
rows <- list()
add <- function(tag, fml, data, family) {
  rows[[length(rows) + 1L]] <<- list(tag = tag, fml = fml, data = data,
                                     family = family)
}
if (set == "sweep122") {
  spec <- switch(part,
    gauss = list(f = "gaussian", extra = "approx"),
    pb = list(f = c("poisson", "bernoulli"), extra = ""),
    student = list(f = "student", extra = "big"))
  for (fam in spec$f) for (n in c(40L, 120L)) for (seed in 1:2) {
    for (body in gf_bodies) {
      add(sprintf("%s_%s_n%d_s%d", body, fam, n, seed),
          gf_sweep_spell(body), gf_sweep_mk(body, fam, n, seed), gf_famf(fam))
    }
  }
  if (spec$extra == "big") {
    for (body in c("exp_a", "exp_b", "pow_a", "mm_K")) {
      add(sprintf("%s_gaussian_n400_s2", body), gf_sweep_spell(body),
          gf_sweep_mk(body, "gaussian", 400L, 2L), gaussian())
    }
  }
  if (spec$extra == "approx") {
    for (body in gf_bodies) {
      add(sprintf("%s_gaussian_n120_s2_k12", body), gf_sweep_spell(body, 12L),
          gf_sweep_mk(body, "gaussian", 120L, 2L), gaussian())
    }
  }
} else if (set == "fresh") {
  source(file.path(GF_WT, "dev/review-fitedge-designs.R"))
  g <- expand.grid(seed = 1:2, n = c(60L, 120L), ci = 1:3,
                   design = rv_designs, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(g))) {
    r <- g[i, ]
    sp <- rv_spell[[r$design]]
    add(sprintf("%s %s %d %d", r$design, c("sin", "wig", "exp")[r$ci], r$n,
                r$seed), sp$nl, rv_mk(r$design, r$seed, r$n, r$ci), sp$fam)
  }
} else if (set == "lane96") {
  g <- expand.grid(seed = 1:3, n = c(60L, 120L), ci = 1:4,
                   design = names(gf_l96_spell), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(g))) {
    r <- g[i, ]
    sp <- gf_l96_spell[[r$design]]
    add(sprintf("%s %s %d %d", r$design, gf_l96_curves[r$ci], r$n, r$seed),
        sp$nl, gf_l96_mk(r$seed, r$n, r$ci, r$design), sp$fam)
  }
} else if (set == "ng150") {
  for (cs in strsplit(part, ",")[[1]]) for (cn in names(gf_ng_curves)) {
    for (n in c(40L, 80L, 160L)) {
      m <- gf_ng_mk(cs, cn, n, 1L)
      add(sprintf("%s %s %d 1", cs, cn, n), m$fml, m$d, m$fam)
    }
  }
}
cat("ROWS", length(rows), "\n")
i <- 0L
for (r in rows) {
  i <- i + 1L
  arms <- if (i %% 2) c("lane", "base") else c("base", "lane")
  for (arm in arms) {
    arm_set(arm)
    t0 <- proc.time()[[3]]
    q <- gf_fit(r$fml, data = r$data, family = r$family)
    wall <- proc.time()[[3]] - t0
    f <- q$fit
    if (inherits(f, "error")) {
      cat(sprintf("ROW %s | %s ERROR %s\n", r$tag, arm,
                  substr(conditionMessage(f), 1, 100)))
      next
    }
    gs <- f$opt$gp_second
    th <- f$opt$par[names(f$opt$par) == "theta"]
    cat(sprintf(paste0("ROW %s | %s ll %.6f conv %d nwarn %d wall %.2f ",
                       "runs %s gain %s theta %s | %s\n"),
                r$tag, arm, as.numeric(logLik(f)), f$opt$convergence,
                length(q$warn), wall,
                if (is.null(gs)) "NA" else format(gs[["runs"]]),
                if (is.null(gs)) "NA" else format(signif(gs[["gain"]], 4)),
                paste(sprintf("%.3f", th), collapse = ","),
                paste(unique(sub(":.*", "", q$warn)), collapse = ";")))
    flush.console()
  }
}
cat("DONE\n")
