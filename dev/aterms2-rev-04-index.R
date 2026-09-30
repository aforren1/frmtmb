# Reviewer, claim 3: index() and mi(x, idx = ) against brms's standata:
# idxl, Jmi and N per response, and the refusals. Seed 404.
# Log: dev/aterms2-rev-log-04-index.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) suppressWarnings(suppressMessages(expr))
err <- function(expr) tryCatch({q(expr); "no error"},
                               error = function(e) conditionMessage(e))
mi_idxl <- function(f) {
  for (lp in f$frame$linpreds) for (mt in lp$mi) {
    if (!is.null(mt$idxl)) return(as.integer(mt$idxl))
  }
  NULL
}
compare <- function(label, fform, bform, data) {
  cat("\n==", label, "\n")
  f <- tryCatch(q(eval(bquote(frm(.(fform), data = data,
                                  family = gaussian())))),
                error = function(e) conditionMessage(e))
  s <- tryCatch(q(brms::standata(bform, data = data)),
                error = function(e) conditionMessage(e))
  if (is.character(f) || is.character(s)) {
    cat("  frmtmb:", if (is.character(f)) f else "fits", "\n")
    cat("  brms  :", if (is.character(s)) s else "standata ok", "\n")
    return(invisible(NULL))
  }
  il <- mi_idxl(f)
  bl <- s[[grep("^idxl", names(s), value = TRUE)[1]]]
  cat("  idxl identical:", identical(il, as.integer(bl)), " n", length(il),
      "\n")
  cat("  N per response frmtmb:", vapply(f$frame$y, NROW, 1L),
      " brms:", unlist(s[grep("^N_", names(s))]), "\n")
  # Jmi: the positions of x's missing values among x's rows
  jm <- f$frame$mi_map$x$rows
  cat("  Jmi identical:", identical(as.integer(jm),
                                     as.integer(s$Jmi_x)), "\n")
  invisible(f)
}
set.seed(404)
n <- 60
d <- data.frame(w = rnorm(n), id = sample(1000, n),
                s = rep(c(TRUE, FALSE, TRUE), length.out = n))
d$x <- rnorm(n)
d$x[c(1, 7, 25)] <- NA
xrows <- which(d$s)
d$ref <- d$id[sample(xrows, n, TRUE)]
d$y <- 1 + 0.5 * d$x[match(d$ref, d$id)] + d$w + rnorm(n)
d$y[is.na(d$y)] <- rnorm(sum(is.na(d$y)))

compare("shuffled numeric ids, x subsetted",
        quote(bf(y ~ mi(x, idx = ref) + w) +
                bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d)
# factor ids whose level orders differ between idx and index
d$idf <- factor(paste0("k", d$id), levels = rev(paste0("k", d$id)))
d$reff <- factor(paste0("k", d$ref))
compare("factor ids, different level orders",
        quote(bf(y ~ mi(x, idx = reff) + w) +
                bf(x | mi() + index(idf) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = reff) + w) +
          brms::bf(x | mi() + index(idf) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d)
# integer index against double idx
d$idd <- as.double(d$id)
compare("integer index, double idx",
        quote(bf(y ~ mi(x, idx = ref) + w) +
                bf(x | mi() + index(idd) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref) + w) +
          brms::bf(x | mi() + index(idd) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d)
# a row dropped for everyone (NA in w) shifts the row numbers
d2 <- d
d2$w[c(4, 30)] <- NA
compare("global NA drop shifts rows",
        quote(bf(y ~ mi(x, idx = ref) + w) +
                bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d2)
# y subsetted too, idx NA outside y's rows
d3 <- d
d3$sy <- rep(c(TRUE, TRUE, FALSE, TRUE), length.out = n)
d3$ref[which(!d3$sy)[1:3]] <- NA
compare("both subsetted, idx NA outside y's rows",
        quote(bf(y | subset(sy) ~ mi(x, idx = ref) + w) +
                bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y | subset(sy) ~ mi(x, idx = ref) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d3)
# duplicated index values among x's rows: both refuse
d4 <- d
d4$id[xrows[2]] <- d4$id[xrows[1]]
compare("duplicated index among x's rows",
        quote(bf(y ~ mi(x, idx = ref) + w) +
                bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d4)
# duplicated only among rows x leaves out: brms checks x's own rows?
d5 <- d
nx <- which(!d5$s)
d5$id[nx[2]] <- d5$id[nx[1]]
compare("duplicated index only outside x's rows",
        quote(bf(y ~ mi(x, idx = ref) + w) +
                bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d5)
# idx value missing from x's index
d6 <- d
d6$ref[3] <- d6$id[nx[1]]
compare("idx value not among x's rows",
        quote(bf(y ~ mi(x, idx = ref) + w) +
                bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d6)
# no index() on x
compare("idx without index()",
        quote(bf(y ~ mi(x, idx = ref) + w) + bf(x | mi() + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref) + w) +
          brms::bf(x | mi() + subset(s) ~ 1) + brms::set_rescor(FALSE), d)
# subsetted x, no idx
compare("subsetted x, mi() without idx",
        quote(bf(y ~ mi(x) + w) + bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d)
# subsetted y, x full, no idx
compare("subsetted y, mi() without idx",
        quote(bf(y | subset(s) ~ mi(x) + w) + bf(x | mi() ~ 1)),
        brms::bf(y | subset(s) ~ mi(x) + w) + brms::bf(x | mi() ~ 1) +
          brms::set_rescor(FALSE), d)
# idx given as an expression
compare("idx = an expression",
        quote(bf(y ~ mi(x, idx = ref + 0) + w) +
                bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref + 0) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d)
# the same x read through two idx variables
d$ref2 <- d$id[sample(xrows, n, TRUE)]
compare("two idx variables for one x",
        quote(bf(y ~ mi(x, idx = ref) + mi(x, idx = ref2) + w) +
                bf(x | mi() + index(id) + subset(s) ~ 1)),
        brms::bf(y ~ mi(x, idx = ref) + mi(x, idx = ref2) + w) +
          brms::bf(x | mi() + index(id) + subset(s) ~ 1) +
          brms::set_rescor(FALSE), d)

cat("\n== predictions on newdata with idx\n")
f <- q(frm(bf(y ~ mi(x, idx = ref) + w) +
             bf(x | mi() + index(id) + subset(s) ~ 1), data = d,
           family = gaussian()))
nd <- d[!is.na(d$x), ]
nd <- nd[sample(nrow(nd)), ]
nd$ref <- nd$id[sample(which(nd$s), nrow(nd), TRUE)]
fy <- fitted(f, newdata = nd, resp = "y")[, "Estimate"]
b <- fixef(f)[, "Estimate"]
man <- b["y_Intercept"] + b["y_mixidxEQref"] * nd$x[match(nd$ref, nd$id)] +
  b["y_w"] * nd$w
cat("  fixef names:", names(b), "\n")
cat("  fitted(newdata) vs manual, max rel:",
    max(abs(fy - man)) / max(abs(man)), "\n")
