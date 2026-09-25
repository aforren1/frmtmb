# Lane wt-mvprior: the base problems the punch-round reviewer found near
# this area, each run through frmtmb (lane library) and brms 2.23.0's
# stancode() on the same data, to file them with a measurement.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressMessages(library(frmtmb))
set.seed(4)
n <- 60
d <- data.frame(x = rnorm(n), w = rnorm(n), xe = rnorm(n), sde = 0.3)
d$y1 <- d$x + rnorm(n)
d$y2 <- d$x + rnorm(n)
d$y3 <- d$w + rnorm(n)
d$o <- factor(cut(d$x + rlogis(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE),
              ordered = TRUE)
res <- function(expr) {
  tryCatch({
    force(expr)
    "ACCEPT"
  }, error = function(e) paste("REFUSE:", gsub("[[:space:]]+", " ",
                                                 conditionMessage(e))))
}
cases <- list(
  three_bf = quote(bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ w) + set_rescor(FALSE)),
  mv_cumulative = quote(bf(y1 ~ x) + bf(o ~ x, family = cumulative()) +
                          set_rescor(FALSE)),
  me = quote(bf(y1 ~ me(xe, sde))),
  zero_intercept = quote(bf(y1 ~ 0 + Intercept + x)),
  student_rescor = quote(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)))
fam <- list(three_bf = gaussian(), mv_cumulative = NULL, me = gaussian(),
            zero_intercept = gaussian(), student_rescor = NULL)
for (nm in names(cases)) {
  fm <- if (nm == "student_rescor") student() else fam[[nm]]
  cat(sprintf("frmtmb %-15s %s\n", nm,
              res(suppressWarnings(frm(eval(cases[[nm]],
                                            asNamespace("frmtmb")),
                                       data = d, family = fm)))))
}
cat(sprintf("frmtmb %-15s %s\n", "mvbf(3)",
            res(frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ w)) +
                      set_rescor(FALSE) + gaussian(), data = d))))
suppressMessages(library(brms))
for (nm in names(cases)) {
  fm <- if (nm == "student_rescor") brms::student() else fam[[nm]]
  cat(sprintf("brms   %-15s %s\n", nm,
              res(brms::stancode(eval(cases[[nm]], asNamespace("brms")),
                                 data = d, family = fm))))
}
