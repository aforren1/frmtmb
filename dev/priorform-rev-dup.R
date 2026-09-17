# Reviewer, lane wt-priorform: the duplicated group-level refusal on the
# designs the worker's generated set did not hold (gr(), known and
# structured covariances, mm(), |ID|, several dpars and responses,
# nonlinear parameters) and on cs() as a covariance structure.
#   Rscript dev/priorform-rev-dup.R brms|ref|lane
# Writes dev/priorform-rev-dup-<mode>.tsv. Seed 20260916.
mode <- commandArgs(trailingOnly = TRUE)[1]
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
pin <- "C:/Users/adf44/source/r/pinlib"
lib <- switch(mode,
  brms = c(pin, user),
  ref = c("C:/Users/adf44/source/r/rellib-r3", pin, user),
  lane = c("C:/Users/adf44/source/r/priorform-lib", pin, user))
.libPaths(lib)
if (mode == "brms") {
  suppressMessages(library(brms))
} else {
  suppressMessages(library(frmtmb))
  cat("frmtmb from", find.package("frmtmb"), "\n")
}
set.seed(20260916)
ng <- 12; per <- 8; n <- ng * per
lev <- paste0("i", 1:ng)
A <- diag(ng); for (i in seq(1, ng - 1, 2)) A[i, i + 1] <- A[i + 1, i] <- 0.5
dimnames(A) <- list(lev, lev)
W <- matrix(0, ng, ng); for (i in 1:(ng - 1)) W[i, i + 1] <- W[i + 1, i] <- 1
dimnames(W) <- list(lev, lev)
V3 <- diag(3) * 0.3; V3[1, 2] <- V3[2, 1] <- 0.1
d <- data.frame(g = factor(rep(lev, each = per), levels = lev),
                h = factor(rep(1:per, ng)),
                x = rnorm(n), z = rnorm(n))
d$id <- d$g
d$id2 <- d$g
d$g1 <- factor(sample(lev, n, TRUE), levels = lev)
d$g2 <- factor(sample(lev, n, TRUE), levels = lev)
d$f <- factor(rep(c("a", "b", "c"), length.out = n))
d$time <- factor(rep(1:4, length.out = n))
if (mode != "brms") d$pos <- num_factor(rep(1:4, length.out = n))
d$y <- rnorm(n) + rnorm(ng)[d$g]
d$y1 <- d$y + rnorm(n); d$y2 <- rnorm(n)
d$yo <- factor(cut(d$y, 4, labels = FALSE), ordered = TRUE)
d$cnt <- rpois(n, 3)
d2 <- list(A = A, W = W, V3 = V3)

# id, family, and the call text; B marks a design brms can express
designs <- list(
  list("gr_cov+plain", "g", "bf(y ~ (1 | gr(id, cov = A)) + (1 | id))", TRUE),
  list("gr_cov+copy", "g", "bf(y ~ (1 | gr(id, cov = A)) + (1 | id2))", TRUE),
  list("gr_cov_slope+plain_slope", "g",
       "bf(y ~ (x | gr(id, cov = A)) + (0 + x | id))", TRUE),
  list("gr_cov+gr_cov", "g",
       "bf(y ~ (1 | gr(id, cov = A)) + (1 | gr(id, cov = A)))", TRUE),
  list("gr_by+plain", "g", "bf(y ~ (1 | gr(g, by = f)) + (1 | g))", TRUE),
  list("gr_student+plain", "g",
       "bf(y ~ (1 | gr(g, dist = 'student')) + (1 | g))", TRUE),
  list("g+g:h", "g", "bf(y ~ (1 | g) + (1 | g:h))", TRUE),
  list("g+g/h", "g", "bf(y ~ (1 | g) + (1 | g/h))", TRUE),
  list("g:h+h:g", "g", "bf(y ~ (1 | g:h) + (1 | h:g))", TRUE),
  list("int+0slope", "g", "bf(y ~ (1 | g) + (0 + x | g))", TRUE),
  list("dblbar", "g", "bf(y ~ (1 + x || g))", TRUE),
  list("dblbar+int", "g", "bf(y ~ (1 + x || g) + (1 | g))", TRUE),
  list("dblbar_f", "g", "bf(y ~ (0 + f || g) + (1 | g))", TRUE),
  list("factor_slope+int", "g", "bf(y ~ (0 + f | g) + (1 | g))", TRUE),
  list("factor+int", "g", "bf(y ~ (f | g) + (1 | g))", TRUE),
  list("mm+member", "g", "bf(y ~ (1 | mm(g1, g2)) + (1 | g1))", TRUE),
  list("mm+mm", "g", "bf(y ~ (1 | mm(g1, g2)) + (1 | mm(g1, g2)))", TRUE),
  list("mm+mm_swapped", "g", "bf(y ~ (1 | mm(g1, g2)) + (1 | mm(g2, g1)))", TRUE),
  list("ID_int+ID_slope0", "g", "bf(y ~ (1 | p | g) + (0 + x | p | g))", TRUE),
  list("ID_int+ID_int", "g", "bf(y ~ (1 | p | g) + (x | p | g))", TRUE),
  list("ID_int+q_int", "g", "bf(y ~ (1 | p | g) + (x | q | g))", TRUE),
  list("mu_sigma_same", "g", "bf(y ~ (1 | g), sigma ~ (1 | g))", TRUE),
  list("mu_sigma_ID", "g", "bf(y ~ (1 | p | g), sigma ~ (1 | p | g))", TRUE),
  list("sigma_dup", "g", "bf(y ~ 1, sigma ~ (1 | g) + (x | g))", TRUE),
  list("mv_same", "mv", "mvbf(bf(y1 ~ (1 | g)), bf(y2 ~ (1 | g)))", TRUE),
  list("mv_ID", "mv", "mvbf(bf(y1 ~ (1 | p | g)), bf(y2 ~ (1 | p | g)))", TRUE),
  list("mvbind", "mv", "bf(mvbind(y1, y2) ~ (1 | g))", TRUE),
  list("mv_dup_one", "mv", "mvbf(bf(y1 ~ (1 | g) + (1 | g)), bf(y2 ~ 1))", TRUE),
  list("nl_same_g", "g",
       "bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1 + (1 | g), nl = TRUE)", TRUE),
  list("nl_dup_in_a", "g",
       "bf(y ~ a + b * x, a ~ 1 + (1 | g) + (1 | g), b ~ 1, nl = TRUE)", TRUE),
  list("nl_gr_cov_a_plain_b", "g",
       "bf(y ~ a + b * x, a ~ 1 + (1 | gr(id, cov = A)), b ~ 1 + (1 | id), nl = TRUE)", TRUE),
  list("car+int", "g", "bf(y ~ car(W, gr = g) + (1 | g))", TRUE),
  list("smooth_re+int", "g", "bf(y ~ s(g, bs = 're') + (1 | g))", TRUE),
  list("ordinal_cs+re", "o", "bf(yo ~ cs(x) + (1 | g))", TRUE),
  list("gp+int", "g", "bf(y ~ gp(x) + (1 | g))", TRUE),
  # frmtmb's glmmTMB structures; brms has none of these
  list("ar1+int", "g", "bf(y ~ ar1(0 + time | g) + (1 | g))", FALSE),
  list("ar1+diag", "g", "bf(y ~ ar1(0 + time | g) + diag(0 + time | g))", FALSE),
  list("us+diag", "g", "bf(y ~ us(0 + time | g) + diag(0 + time | g))", FALSE),
  list("rr+diag", "g", "bf(y ~ rr(0 + time | g, d = 1) + diag(0 + time | g))", FALSE),
  list("cs_cov+int", "g", "bf(y ~ cs(0 + time | g) + (1 | g))", FALSE),
  list("cs_cov+diag", "g", "bf(y ~ cs(0 + time | g) + diag(0 + time | g))", FALSE),
  list("homcs+diag", "g", "bf(y ~ homcs(0 + time | g) + diag(0 + time | g))", FALSE),
  list("equalto+us", "g", "bf(y ~ equalto(0 + f | g, V3) + (0 + f | g))", FALSE),
  list("equalto+diag", "g", "bf(y ~ equalto(0 + f | g, V3) + diag(0 + f | g))", FALSE),
  list("exp+diag", "g", "bf(y ~ exp(0 + pos | g) + diag(0 + pos | g))", FALSE),
  list("ou+int", "g", "bf(y ~ ou(0 + pos | g) + (1 | g))", FALSE),
  list("gr_prec+plain", "g", "bf(y ~ (1 | gr(id, prec = A)) + (1 | id))", FALSE),
  list("propto+us", "g", "bf(y ~ propto(0 + f | g, V3) + (0 + f | g))", FALSE),
  # cs() as a covariance structure against the category-specific branch
  list("cs_cov_ordinal", "o", "bf(yo ~ x + cs(0 + time | g))", FALSE),
  list("cs_cov_dblbar", "g", "bf(y ~ cs(0 + time || g))", FALSE),
  list("cs_cov_twice", "g", "bf(y ~ cs(0 + time | g) + cs(0 + time | h))", FALSE),
  list("cs_catspec_and_cov", "o", "bf(yo ~ cs(x) + cs(0 + time | g))", FALSE),
  list("cs_catspec_selfint", "o", "bf(yo ~ cs(x) * cs(x))", TRUE),
  list("cs_catspec_cov_cross", "o", "bf(yo ~ cs(x) * cs(0 + time | g))", FALSE)
)

fam_of <- function(k) switch(k, g = "gaussian()", o = "sratio()",
                             mv = NULL)
one <- function(ds) {
  id <- ds[[1]]; fk <- ds[[2]]; txt <- ds[[3]]; brms_ok <- ds[[4]]
  if (mode == "brms" && !brms_ok) return(c(id, "n/a", ""))
  res <- tryCatch({
    if (mode == "brms") {
      form <- eval(parse(text = txt))
      fam <- if (is.null(fam_of(fk))) NULL else
        eval(parse(text = fam_of(fk)))
      if (is.null(fam)) {
        suppressMessages(stancode(form, data = d, data2 = d2))
      } else {
        suppressMessages(stancode(form, data = d, family = fam, data2 = d2))
      }
    } else {
      form <- eval(parse(text = txt))
      fam <- if (is.null(fam_of(fk))) NULL else
        eval(parse(text = fam_of(fk)))
      if (!is.null(fam)) form <- form + fam
      suppressWarnings(suppressMessages(
        frm(form, data = d, data2 = d2, dry_run = "frame")))
    }
    c("accept", "")
  }, error = function(e) {
    msg <- gsub("[\t\r\n]+", " ", conditionMessage(e))
    v <- if (grepl("Duplicated group-level", msg)) "refuse-dup" else
      if (grepl("is invalid", msg)) "refuse-term" else "error-other"
    c(v, substr(msg, 1, 160))
  })
  c(id, res)
}
out <- do.call(rbind, lapply(designs, one))
colnames(out) <- c("id", "verdict", "message")
path <- sprintf("C:/Users/adf44/source/r/frmtmb-wt-priorform/dev/priorform-rev-dup-%s.tsv", mode)
write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
print(table(out[, "verdict"]))
