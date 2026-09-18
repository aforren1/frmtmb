.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")

set.seed(20260917)
d <- data.frame(g1 = rep(c(1, 1, 2, 2), each = 6),
                g2 = rep(c(1, 2, 1, 2), each = 6),
                t = c(1:6, 1:6, 1:6, 7:12))
d$x <- rnorm(24)
d$g <- interaction(d$g1, d$g2)
d$y <- d$x + as.numeric(stats::filter(rnorm(24), 0.6, "recursive"))
d$gf1 <- factor(d$g1)
d$gf2 <- factor(d$g2)
d$g3 <- rep(c(1, 2), 12)
d$g.dot <- d$g1
d[[".lead"]] <- d$g1
d[["my var"]] <- d$g1
d$t.dot <- d$t
d[["my t"]] <- d$t
gv <- d$g1
gfac <- factor(d$g1)

cases <- c(
  "ar(t, g1, cov = TRUE)",
  "ar(t, factor(g1), cov = TRUE)",
  "ar(t, interaction(g1, g2), cov = TRUE)",
  "ar(t, g1:g2, cov = TRUE)",
  "ar(t, gf1:gf2, cov = TRUE)",
  "ar(t, g1:g2:g3, cov = TRUE)",
  "ar(t, g.dot, cov = TRUE)",
  "ar(t, .lead, cov = TRUE)",
  "ar(t, `my var`, cov = TRUE)",
  "ar(t, d$g1, cov = TRUE)",
  "ar(t, gv, cov = TRUE)",
  "ar(t, gfac, cov = TRUE)",
  "ar(t, get(\"g1\"), cov = TRUE)",
  "ar(t, g1/g2, cov = TRUE)",
  "ar(t, g1 + g2, cov = TRUE)",
  "ar(t, g1 * g2, cov = TRUE)",
  "ar(t, (g1), cov = TRUE)",
  "ar(t, 1, cov = TRUE)",
  "ar(t, \"g1\", cov = TRUE)",
  "ar(t, NA, cov = TRUE)",
  "ar(x + t, g, cov = TRUE)",
  "ar(t.dot, g, cov = TRUE)",
  "ar(`my t`, g, cov = TRUE)",
  "ar(d$t, g, cov = TRUE)",
  "ar(as.integer(t), g, cov = TRUE)",
  "ar(1, g, cov = TRUE)",
  "ar(\"t\", g, cov = TRUE)",
  "ar(NA, g, cov = TRUE)",
  "ar(gr = g, cov = TRUE)",
  "cosy(gr = g1:g2)",
  "cosy(gr = factor(g1))",
  "unstr(t, g1:g2)",
  "ma(t, g1:g2, cov = TRUE)",
  "arma(t, g1:g2, cov = TRUE)"
)

short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)), 1, 110)

cat("CASE\tBRMS_CTOR\tBRMS_STANDATA\n")
for (cs in cases) {
  v1 <- tryCatch({
    eval(parse(text = paste0("brms::", cs)))
    "OK"
  }, error = function(e) paste0("REFUSE: ", short(e)))
  v2 <- tryCatch({
    sd <- suppressWarnings(suppressMessages(
      brms::make_standata(stats::as.formula(paste("y ~ x +", cs)), d,
                          family = brms::brmsfamily("gaussian"))))
    paste0("OK n=", sd$N, " nJ=",
           if (!is.null(sd$J_lag)) length(unique(sd$J_lag)) else NA)
  }, error = function(e) paste0("REFUSE: ", short(e)))
  cat(cs, "\t", v1, "\t", v2, "\n", sep = "")
}
