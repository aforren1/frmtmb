# Review of lane prilabel: brms 2.23.0 against the lane's frmtmb on the
# cases the lane did not list. No compile. Data dev/dblbar-data.R seed
# 20261013; extra columns seed 7 (as dev/prilabel-sdcmp.R) and seed 8.
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
})
dir.create(file.path(wt, "dev/prilabel-rev-log"), showWarnings = FALSE)
sink(file.path(wt, "dev/prilabel-rev-log/probe.txt"), split = TRUE)
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "; brms",
    format(packageVersion("brms")), "\n")
dd <- dblbar_data()
set.seed(7)
dd$x1 <- rnorm(nrow(dd))
dd$x2 <- rnorm(nrow(dd))
dd$w1 <- runif(nrow(dd))
dd$w2 <- 1 - dd$w1
dd$g3 <- factor(sample(c("r1", "r2"), nrow(dd), TRUE))
dd$yz <- dd$y + rnorm(nrow(dd))
set.seed(8)
# a by-variable constant within g1
dd$g1b <- factor(ifelse(as.integer(dd$g1) <= 3L, "lo", "hi"))

key_cols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
brms_rows <- function(bf_b, data = dd) {
  p <- as.data.frame(brms::default_prior(bf_b, data = data))
  p <- p[p$class %in% c("sd"), key_cols]
  rownames(p) <- NULL
  p
}
frm_rows <- function(bf_f) {
  p <- as.data.frame(frmtmb::default_prior(bf_f, data = dd))
  p <- p[p$class %in% c("sd"), key_cols]
  rownames(p) <- NULL
  p
}
try_ <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
# brms: the Stan lprior lines for sd_ and the frame_re map
brms_sd_lines <- function(bf_b, prior) {
  sc <- brms::stancode(bf_b, data = dd, prior = prior)
  grep("lprior.*sd_", strsplit(sc, "\n")[[1]], value = TRUE)
}
brms_re <- function(bf_b) {
  re <- as.data.frame(brms:::frame_re(brms:::brmsterms(bf_b), data = dd))
  re[, intersect(c("id", "group", "gn", "coef", "resp", "dpar", "nlpar",
                   "cor", "by"), names(re))]
}
# frmtmb: the entries a prior resolves to, with the block and column
frm_entries <- function(bf_f, prior) {
  des <- frmtmb:::prior_design(bf_f, dd, NULL, list(),
                               check_response = FALSE)
  ent <- frmtmb:::resolve_priorlist(des, prior)$entries
  out <- NULL
  for (e in ent) {
    if (!identical(e$comp, "theta")) next
    for (bk in des$frame$re_blocks) {
      j <- match(e$idx, bk$theta_idx)
      if (is.na(j)) next
      out <- rbind(out, data.frame(
        theta = e$idx, group = bk$group_name, cs = bk$covstruct,
        pos = j, kind = e$dist$kind,
        scale = if (is.null(e$dist$scale)) NA else e$dist$scale,
        cnms = paste(bk$cnms, collapse = ",")))
    }
  }
  out
}
mk <- function(pkg, f, s = NULL) {
  fn <- getExportedValue(pkg, "bf")
  if (is.null(s)) fn(f) else fn(f, s)
}
cmp_rows <- function(f, s = NULL) {
  lab <- paste(deparse1(f), if (!is.null(s)) deparse1(s))
  b <- try_(brms_rows(mk("brms", f, s)))
  r <- try_(frm_rows(mk("frmtmb", f, s)))
  cat(sprintf("\nROWS %s %s\n", if (identical(b, r)) "SAME" else "DIFF",
              lab))
  if (!identical(b, r)) {
    cat(" brms:\n"); print(b)
    cat(" frmtmb:\n"); print(r)
  } else {
    print(b)
  }
}

cat("\n######## 1. default_prior sd rows on unlisted cases\n")
cmp_rows(y ~ x + (1 + x | g1/g2/g3))
cmp_rows(y ~ x + (1 + x | g1) + (1 + x | g1:g2))
cmp_rows(y ~ x + (1 + x | mm(m1, m2, weights = cbind(w1, w2))))
cmp_rows(y ~ x + (1 + x | gr(g1, by = g1b)))
cmp_rows(y ~ x + (1 + f || g))
cmp_rows(y ~ x + (1 | g1:g2/g3))
cmp_rows(y ~ x + (1 | g1/g2), sigma ~ (1 | g1/g2))
cmp_rows(y ~ x + (1 + x | g2/g1))

cat("\n######## 2. mvbf with |ID| across responses\n")
mvb <- brms::bf(y ~ x + (1 + x | p | g)) + brms::bf(yz ~ z + (1 | p | g)) +
  brms::set_rescor(FALSE)
mvf <- frmtmb::bf(y ~ x + (1 + x | p | g)) + frmtmb::bf(yz ~ z + (1 | p | g)) +
  frmtmb::set_rescor(FALSE)
b <- try_(brms_rows(mvb))
r <- try_(frm_rows(mvf))
cat("ROWS", if (identical(b, r)) "SAME" else "DIFF", "\n")
print(b)
if (!identical(b, r)) print(r)
for (a in list(list(group = "g", coef = "x", resp = "y"),
               list(group = "g", coef = "Intercept", resp = "yz"),
               list(group = "g", resp = "yz"))) {
  cat("\n-- prior", paste(names(a), unlist(a), sep = "=", collapse = ", "),
      "\n")
  pb <- do.call(brms::set_prior, c(list("normal(0, 0.5)", class = "sd"), a))
  print(try_(brms_sd_lines(mvb, pb)))
  pf <- do.call(frmtmb::set_prior, c(list("normal(0, 0.5)", class = "sd"), a))
  print(try_(frm_entries(mvf, pf)))
}
print(brms_re(mvb))

cat("\n######## 3. per-coef priors on the unlisted cases\n")
pcs <- list(
  list(y ~ x + (1 + x | g1/g2/g3), list(group = "g1:g2:g3", coef = "x")),
  list(y ~ x + (1 + x | g1/g2/g3), list(group = "g1:g2", coef = "Intercept")),
  list(y ~ x + (1 + x | mm(m1, m2, weights = cbind(w1, w2))),
       list(group = "mmm1m2", coef = "x")),
  list(y ~ x + (1 + x | gr(g1, by = g1b)), list(group = "g1", coef = "x")),
  list(y ~ x + (1 + x | gr(g1, by = g1b)),
       list(group = "g1", coef = "Intercept")),
  list(y ~ x + (1 + x | gr(g1, by = g1b)), list(group = "g1")),
  list(y ~ x + (1 + f || g), list(group = "g", coef = "fc")),
  list(y ~ x + (1 + x | g1) + (1 + x | g1:g2), list(group = "g1:g2",
                                                     coef = "x"))
)
for (pc in pcs) {
  a <- pc[[2]]
  cat("\n--", deparse1(pc[[1]]), "::",
      paste(names(a), unlist(a), sep = "=", collapse = ", "), "\n")
  pb <- do.call(brms::set_prior, c(list("normal(0, 0.5)", class = "sd"), a))
  print(try_(brms_sd_lines(brms::bf(pc[[1]]), pb)))
  pf <- do.call(frmtmb::set_prior, c(list("normal(0, 0.5)", class = "sd"), a))
  print(try_(frm_entries(frmtmb::bf(pc[[1]]), pf)))
}
print(brms_re(brms::bf(y ~ x + (1 + x | gr(g1, by = g1b)))))

cat("\n######## 4. duplicates and odd nests\n")
for (f in list(y ~ x + (1 | g1/g2) + (1 | g1),
               y ~ x + (1 | g1/g2) + (0 + x | g1),
               y ~ x + (1 | g1/g2) + (1 | g1:g2),
               y ~ x + (1 | p | g1/g2),
               y ~ x + (1 | g1/g2) + (1 | g2:g1),
               y ~ x + (1 | g1:g2) + (1 | g2:g1))) {
  cat("\n--", deparse1(f), "\n")
  cat(" brms:  ", try_(paste(unique(brms_rows(brms::bf(f))$group),
                             collapse = " ; ")), "\n")
  cat(" frmtmb:", try_(paste(unique(frm_rows(frmtmb::bf(f))$group),
                             collapse = " ; ")), "\n")
}

cat("\n######## 5. refusals where the condition is absent\n")
pr <- frmtmb::set_prior("normal(0, 0.5)", class = "sd", group = "g",
                        coef = "x")
for (f in list(
  list(y ~ x + (1 + x | g), sigma ~ homdiag(1 + x | g)),
  list(y ~ x + homdiag(1 + x | g), sigma ~ (1 + x | g)),
  list(y ~ x + homdiag(1 + z | g) + (0 + x | g)),
  list(y ~ x + (1 + x | g) + homdiag(1 + x | h)),
  list(y ~ x + ar1(0 + f | g) + (0 + x | g)))) {
  cat("\n--", paste(vapply(f, deparse1, ""), collapse = " ; "), "\n")
  bff <- do.call(frmtmb::bf, f)
  print(try_(frm_entries(bff, pr)))
  if (length(f) == 2L) {
    print(try_(frm_entries(bff, frmtmb::set_prior(
      "normal(0, 0.5)", class = "sd", group = "g", coef = "x",
      dpar = "sigma"))))
  }
}
cat("\n-- brms coef on (1 + x | g) with sigma ~ (1 + x | g)\n")
print(brms_sd_lines(brms::bf(y ~ x + (1 + x | g), sigma ~ (1 + x | g)),
                    brms::set_prior("normal(0, 0.5)", class = "sd",
                                    group = "g", coef = "x")))
cat("\n-- a Student-t block per coef, and a quadrature fit of one t term\n")
print(try_(frm_entries(frmtmb::bf(y ~ x + (1 + x | gr(g, dist = "student"))),
                       pr)))
print(brms_sd_lines(brms::bf(y ~ x + (1 + x | gr(g, dist = "student"))),
                    brms::set_prior("normal(0, 0.5)", class = "sd",
                                    group = "g", coef = "x")))
q1 <- try_(frm(frmtmb::bf(y ~ x + (1 | gr(g, dist = "student")) + (1 | g1)),
               data = dd, quadrature = TRUE))
cat("quadrature, one t term + gaussian term on another factor:",
    if (is.character(q1)) q1 else "fits", "\n")
q2 <- try_(frm(frmtmb::bf(y ~ x + (1 | gr(g, dist = "student")) +
                            (1 | gr(g1, dist = "student"))),
               data = dd, quadrature = TRUE))
cat("quadrature, t terms on two factors:",
    if (is.character(q2)) q2 else "fits", "\n")
q3 <- try_(frm(frmtmb::bf(y ~ x + (1 | gr(g, dist = "student")) + (0 + x | g)),
               data = dd, quadrature = TRUE))
cat("quadrature, t term + gaussian term on g:",
    if (is.character(q3)) q3 else "fits", "\n")

cat("\n######## 6. adjusted expectations, brms side\n")
data("epilepsy", package = "brms")
vp <- brms::validate_prior(brms::prior(cauchy(0, 2), class = sd,
                                       group = patient),
                           count ~ zAge + (1 | patient), data = epilepsy,
                           family = poisson())
print(as.data.frame(vp)[vp$class == "sd", c("prior", "class", "coef",
                                             "group", "source")])
dB <- data.frame(yb = plogis(rnorm(40)), x = rnorm(40),
                 g = factor(rep(1:4, 10)))
vp2 <- brms::validate_prior(brms::prior(normal(0, 5), class = sd, group = g),
                            brms::bf(yb ~ x + (1 | g), phi ~ (1 | g)),
                            family = brms::Beta(), data = dB)
print(as.data.frame(vp2)[vp2$class == "sd", c("prior", "class", "coef",
                                               "group", "dpar", "source")])
dS <- data.frame(y = rnorm(60), Subject = factor(rep(1:10, 6)),
                 a = factor(rep(1:2, 30)))
print(brms_rows(brms::bf(y ~ 1 + (1 | Subject/a)), data = dS))
dM <- data.frame(y = rnorm(60), x = rnorm(60),
                 g1 = factor(sample(1:6, 60, TRUE)),
                 g2 = factor(sample(1:6, 60, TRUE)))
print(unique(as.data.frame(brms::default_prior(
  y ~ x + (1 | mm(g1, g2)), data = dM))$group))
sink()
