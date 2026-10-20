# Lane prilabel: the frmtmb side of dev/prilabel-brms-ref.R and
# dev/prilabel-brms-fit.R. Arg 1: "base" (rellib-r9) or "lane" (the
# lane library first). Data: dev/dblbar-data.R, seed 20261013.
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
lib <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-prilabel-lib",
         "C:/Users/adf44/source/r/rellib-r9",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
dd <- dblbar_data()
sink(file.path(wt, paste0("dev/prilabel-log/frm-probe-", arm, ".txt")),
     split = TRUE)
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
show_rows <- function(f) {
  cat("\n== default_prior:", deparse1(f$formula %||% f[[1]]), "\n")
  p <- tryCatch(as.data.frame(default_prior(f, data = dd,
                                            family = gaussian())),
                error = function(e) {
                  cat("ERROR:", conditionMessage(e), "\n"); NULL
                })
  if (!is.null(p)) {
    print(p[p$class %in% c("sd", "cor"),
            intersect(c("prior", "class", "coef", "group", "resp", "dpar",
                        "nlpar"), names(p))], row.names = FALSE)
  }
}
cases <- list(
  bf(y ~ x + (1 + x | g)),
  bf(y ~ x + (1 + x || g)),
  bf(y ~ x + (1 + x | gr(g, cor = FALSE))),
  bf(y ~ x + (1 + x | g1/g2)),
  bf(y ~ x + (1 + x | g1:g2)),
  bf(y ~ x + (1 | g1) + (1 | g1:g2)),
  bf(y ~ x + (1 + x | mm(m1, m2))),
  bf(y ~ x + (1 | g), sigma ~ (1 + x | g)),
  bf(y ~ x + (1 | p | g), sigma ~ (1 | p | g)),
  bf(y ~ x + (f | g))
)
for (f in cases) show_rows(f)

try_fit <- function(form, prior) {
  tryCatch({
    fit <- frm(form, family = gaussian(), data = dd, prior = prior)
    pe <- fit$prior_entries %||% NULL
    "accepted"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
}
cat("\n== set_prior(class = sd, group = g, coef = x)\n")
cat(try_fit(bf(y ~ x + (1 + x | g)),
            set_prior("normal(0, 0.5)", class = "sd", group = "g",
                      coef = "x")), "\n")
cat("\n== set_prior(class = sd, group = g1:g2)\n")
cat(try_fit(bf(y ~ x + (1 + x | g1/g2)),
            set_prior("normal(0, 0.5)", class = "sd", group = "g1:g2")),
    "\n")
cat("\n== set_prior(class = sd, group = mmm1m2)\n")
cat(try_fit(bf(y ~ x + (1 + x | mm(m1, m2))),
            set_prior("normal(0, 0.5)", class = "sd", group = "mmm1m2")),
    "\n")

cat("\n== the brms fit's model\n")
fit <- frm(bf(y ~ x + (1 + x | g1/g2) + (1 | mm(m1, m2)), sigma ~ (1 | g)),
           family = gaussian(), data = dd)
v <- variables(fit)
print(grep("^sd_|^cor_", v, value = TRUE))
r <- grep("^r_", v, value = TRUE)
print(r[!duplicated(sub("\\[.*", "", r))])
cat("\n== names(ranef), dimnames heads\n")
re <- ranef(fit)
print(names(re))
for (nm in names(re)) print(head(dimnames(re[[nm]])[[1]], 3))
cat("\n== names(coef)\n")
print(names(coef(fit)))
cat("\n== names(VarCorr)\n")
print(names(VarCorr(fit)))
cat("\n== ngrps\n")
print(tryCatch(ngrps(fit), error = function(e) conditionMessage(e)))
cat("\n== summary random part\n")
s <- summary(fit)
print(names(s$random))
cat("\n== predict, re_formula\n")
nd <- dd[1:3, ]
for (rf in list(~ (1 + x | g1:g2), ~ (1 + x | g1/g2), ~ (1 + x | g2:g1))) {
  p <- tryCatch(dim(fitted(fit, newdata = nd, re_formula = rf)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(deparse1(rf), ":", p, "\n")
}
cat("\n== newdata with a new level of g2\n")
nd2 <- dd[1:2, ]
nd2$g2 <- factor("qNEW")
p4 <- tryCatch(dim(fitted(fit, newdata = nd2, allow_new_levels = TRUE)),
               error = function(e) paste("ERROR:", conditionMessage(e)))
print(p4)
p5 <- tryCatch(fitted(fit, newdata = nd2),
               error = function(e) paste("ERROR:", conditionMessage(e)))
print(p5)
sink()
