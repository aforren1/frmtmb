# Attacks 2 and 6: which block structures the nested route accepts and
# refuses, from Z, and whether each refusal is true, classed, and early.
#   Rscript dev/review-nestquad-struct.R
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
set.seed(611)
ns_ <- 10
d <- do.call(rbind, lapply(seq_len(ns_), function(s) {
  nc <- c(1, 2, 3, 4, 10, 1, 2, 3, 5, 6)[s]
  do.call(rbind, lapply(seq_len(nc), function(c) {
    nr <- sample(1:4, 1)
    data.frame(school = s, cls = c, x = rnorm(nr))
  }))
}))
d$school <- factor(d$school)
# class labels repeat across schools ("1", "2", ...): crossed in names
d$cls_rep <- factor(d$cls)
# unique class labels: nested in Z although written as a crossed sum
d$cls_uni <- factor(paste0(d$school, "_", d$cls))
ua <- rnorm(ns_, 0, 0.8); ui <- rnorm(nlevels(d$cls_uni), 0, 0.8)
d$y <- rbinom(nrow(d), 1, plogis(0.3 + 0.5 * d$x + ua[d$school] +
                                   ui[d$cls_uni]))
d$yg <- 0.3 + 0.5 * d$x + ua[d$school] + ui[d$cls_uni] + rnorm(nrow(d))
d$cnt <- rpois(nrow(d), exp(0.3 + 0.5 * d$x + ua[d$school] * 0.5))
cat("rows", nrow(d), "schools", ns_, "classes", nlevels(d$cls_uni), "\n")

run <- function(tag, fo, data = d, ...) {
  t0 <- proc.time()[[3]]
  r <- try_msg(frm(fo, data = data, quadrature = TRUE, ...))
  el <- proc.time()[[3]] - t0
  if (inherits(r, "error")) {
    cat(sprintf("REFUSE %-22s %5.2fs class=%s\n  msg: %s\n", tag, el,
                paste(class(r), collapse = "/"), conditionMessage(r)))
  } else {
    cat(sprintf("FIT    %-22s %5.2fs logLik %.8f nblocks %d random %d\n",
                tag, el, as.numeric(logLik(r)), length(r$frame$re_blocks),
                length(r$obj$env$random)))
  }
  if (length(attr(r, "warns"))) cat("  warns:", attr(r, "warns"), sep = "\n   ")
  invisible(r)
}
bern <- bernoulli()
# 1. nested via the slash, the reference spelling
f1 <- run("slash", bf(y ~ x + (1 | school/cls_uni)) + bern)
# 2. (1 | school) + (1 | school:cls_rep): same model, other spelling
f2 <- run("colon", bf(y ~ x + (1 | school) + (1 | school:cls_rep)) + bern)
# 3. unique labels written as a crossed sum: Z nested, must fit
f3 <- run("sum-unique", bf(y ~ x + (1 | school) + (1 | cls_uni)) + bern)
# 4. labels repeat across schools, written as a sum: Z crossed, refuse
f4 <- run("sum-repeated", bf(y ~ x + (1 | school) + (1 | cls_rep)) + bern)
# 5. reversed order of the blocks
f5 <- run("reversed", bf(y ~ x + (1 | cls_uni) + (1 | school)) + bern)
for (f in list(f2, f3, f5)) {
  if (!inherits(f, "error") && !inherits(f1, "error")) {
    cat(sprintf("  diff vs slash: %.3e\n",
                as.numeric(logLik(f)) - as.numeric(logLik(f1))))
  }
}
# 6. partially crossed: three classes span two schools
dp <- d
rows <- which(dp$cls_uni %in% c("5_1", "5_2", "5_3"))
half <- rows[seq(1, length(rows), by = 2)]
dp$school[half] <- "4"
cat("partial: rows moved", length(half), "\n")
f6 <- run("partial-cross", bf(y ~ x + (1 | school/cls_uni)) + bern, dp)
f6b <- run("partial-cross-sum", bf(y ~ x + (1 | school) + (1 | cls_uni)) +
             bern, dp)
# 7. a level of one class only per school everywhere: both nest in each
# other
d1 <- d[d$cls == 1, ]
d1$cls_uni <- droplevels(d1$cls_uni)
f7 <- run("one-class-per-school", bf(y ~ x + (1 | school/cls_uni)) + bern, d1)
# 8. two blocks on one factor, same predictor: (1 | g) + (0 + x | g) is
# admitted? (1 + x || g) is one diag block of dim 2 and refused
f8 <- run("int+slope two blocks", bf(y ~ x + (1 | school) + (0 + x | school)) +
            bern)
f8b <- run("int+slope double bar", bf(y ~ x + (1 + x || school)) + bern)
# 9. slope in each block
f9 <- run("slope nested", bf(y ~ x + (0 + x | school/cls_uni)) + bern)
f9b <- run("1+x nested", bf(y ~ x + (1 + x | school/cls_uni)) + bern)
# 10. refusals by name
dd <- d; dd$gc <- factor(paste0(dd$cls_uni, "_", sample(1:2, nrow(dd), TRUE)))
run("three blocks", bf(y ~ x + (1 | school/cls_uni/gc)) + bern, dd)
run("mu+sigma", bf(yg ~ x + (1 | school), sigma ~ (1 | school)) +
      gaussian())
run("mu+sigma nested", bf(yg ~ x + (1 | school), sigma ~ (1 | school:cls_uni)) +
      gaussian())
run("student block", bf(y ~ x + (1 | gr(school, dist = "student")) +
                          (1 | school:cls_uni)) + bern)
run("mv two responses", mvbf(bf(y ~ x + (1 | school/cls_uni)) + bern,
                             bf(yg ~ x + (1 | school/cls_uni)) + gaussian()))
run("mv shared one block", mvbf(bf(y ~ x + (1 | p | school)) + bern,
                             bf(yg ~ x + (1 | p | school)) + gaussian()))
run("nl body", bf(yg ~ a + b * x, a ~ 1 + (1 | school/cls_uni), b ~ 1,
                  nl = TRUE) + gaussian())
d$yo <- factor(cut(d$yg, c(-Inf, -0.5, 0.5, Inf), labels = FALSE),
               ordered = TRUE)
run("cs() acat", bf(yo ~ cs(x) + (1 | school/cls_uni)) + acat())
run("cumulative ok", bf(yo ~ x + (1 | school/cls_uni)) + cumulative())
d$mm1 <- d$school
d$mm2 <- factor(sample(levels(d$school), nrow(d), TRUE),
                levels = levels(d$school))
run("mm + nested", bf(y ~ x + (1 | mm(mm1, mm2)) + (1 | cls_uni)) + bern)
run("grouped mixture", bf(yg ~ x + (1 | school/cls_uni)) +
      mixture(gaussian(), gaussian(), groups = ~school))
# a correlated covariance across levels: gr(cov = ) with a dense matrix
A <- 0.5 * diag(ns_) + 0.5
dimnames(A) <- list(levels(d$school), levels(d$school))
run("gr cov dense", bf(y ~ x + (1 | gr(school, cov = A)) + (1 | cls_uni)) +
      bern, data2 = list(A = A))
# offset and weights pass through
d$off <- runif(nrow(d), -0.3, 0.3); d$w <- runif(nrow(d), 0.5, 1.5)
run("offset+weights", bf(y | weights(w) ~ x + offset(off) +
                           (1 | school/cls_uni)) + bern)
# nested and a by-split
d$by <- factor(c("a", "b")[(as.integer(d$school) %% 2) + 1])
run("by-split outer", bf(y ~ x + (1 | gr(school, by = by)) + (1 | cls_uni)) +
      bern)
cat("DONE\n")
