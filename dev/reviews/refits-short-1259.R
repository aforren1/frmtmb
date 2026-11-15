# Short check of the refits fix round, 1259: character ordinal responses
# against brms 2.23.0 (standata only). Seed 2501, n = 60.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
set.seed(2501)
d1 <- data.frame(x = rnorm(60))
d1$yc <- sample(c("none", "mild", "moderate", "severe"), 60, TRUE)
d1$yd <- sample(c("1", "2", "3", "4"), 60, TRUE)
d1$yo <- factor(d1$yc, levels = c("none", "mild", "moderate", "severe"),
                ordered = TRUE)
d1$yi <- as.integer(d1$yd)
d1$yz <- d1$yi - 1L
fams <- list(cumulative = function() brms::cumulative(),
             sratio = function() brms::sratio(),
             cratio = function() brms::cratio(),
             acat = function() brms::acat())
ffams <- list(cumulative = frmtmb::cumulative, sratio = frmtmb::sratio,
              cratio = frmtmb::cratio, acat = frmtmb::acat)
msg <- function(expr) tryCatch({ expr; "FITTED" },
                               error = function(e) conditionMessage(e))
first <- function(m) sub("^(.*?responses\\.).*$", "\\1", m)
for (fn in names(fams)) {
  for (v in c("yc", "yd")) {
    fo <- as.formula(paste(v, "~ x"))
    b <- msg(brms::standata(fo, data = d1, family = fams[[fn]]()))
    f <- msg(frm(fo, data = d1, family = ffams[[fn]]()))
    cat(sprintf("\n[%s %s]\n brms:   %s\n frmtmb: %s\n first sentence equal: %s\n",
                fn, v, b, f, identical(first(b), first(f))))
  }
  for (v in c("yo", "yi")) {
    fo <- as.formula(paste(v, "~ x"))
    r <- tryCatch(round(fixef(frm(fo, data = d1, family = ffams[[fn]]()),
                              flatten = TRUE), 6),
                  error = function(e) conditionMessage(e))
    cat(sprintf("[%s %s] ", fn, v)); print(r)
  }
}
if (exists("hurdle_cumulative", asNamespace("frmtmb"))) {
  b <- msg(brms::standata(yc ~ x, data = d1, family = brms::hurdle_cumulative()))
  f <- msg(frm(yc ~ x, data = d1, family = frmtmb::hurdle_cumulative()))
  cat(sprintf("\n[hurdle_cumulative yc]\n brms:   %s\n frmtmb: %s\n equal: %s\n",
              b, f, identical(first(b), first(f))))
  b <- msg(brms::standata(yz ~ x, data = d1, family = brms::hurdle_cumulative()))
  cat("[hurdle_cumulative yz integer] brms:", b, "\n")
}
go("NA in character response", {
  d2 <- d1; d2$yd[1:2] <- NA
  frm(yd ~ x, data = d2, family = frmtmb::cumulative()) })
go("all-NA-but-one character", {
  d2 <- d1; d2$yd[] <- NA; d2$yd[1] <- "1"
  frm(yd ~ x, data = d2, family = frmtmb::cumulative()) })
go("character categorical unaffected",
   class(frm(yc ~ x, data = d1, family = frmtmb::categorical())))
go("mvbf ordinal + gaussian with character ordinal", {
  d1$g2 <- rnorm(60)
  frm(bf(yd ~ x) + bf(g2 ~ x) , data = d1,
      family = list(frmtmb::cumulative(), gaussian())) })
