# Review of lane nlstart, short check of the fix round: order independence
# with more than four curved coefficients (the probe then varies the four
# with the largest derivative RMS). One model, one data set, three
# spellings: terms in order, terms reversed (declarations reversed too), and
# nlpars renamed so their alphabetical order differs. Also a body whose
# terms read one covariate, so that keys tie.
# Usage: REV_ARM=after Rscript dev/review-nlstart-fix-order5.R
source("dev/review-nlstart-common.R")
F <- list(gaussian = gaussian(), poisson = poisson(), bernoulli = bernoulli(),
          cumulative = cumulative())
terms6 <- c("p1 * (1 - exp(-p2 * x))", "p3 * (1 - exp(-p4 * z))",
            "p5 * exp(p6 * w)")
terms5 <- c("p1 * (1 - exp(-p2 * x))", "p3 * (1 - exp(-p4 * z))",
            "exp(p5 * w)")
termsT <- c("p1 * exp(p2 * x)", "p3 * (1 - exp(-p4 * x))",
            "p5 * exp(p6 * x)")
spell <- function(terms, np, map, rev = FALSE) {
  tt_ <- if (rev) rev(terms) else terms
  body <- paste(tt_, collapse = " + ")
  for (i in seq_len(np)) body <- gsub(paste0("p", i, "\\b"), map[i], body)
  nms <- map[seq_len(np)]
  if (rev) nms <- rev(nms)
  list(body = body, decl = paste(nms, collapse = " + "))
}
maps <- list(b = paste0("b", 1:6), rename = c("zf", "ye", "xd", "wc", "vb", "ua"))
gen <- function(fam, seed, tie = FALSE) {
  set.seed(seed)
  n <- 400
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3))
  lat <- if (tie) 1.2 * exp(0.4 * d$x) + 1.0 * (1 - exp(-0.7 * d$x)) +
    0.3 * exp(-0.5 * d$x) else
    1.5 * exp(0.6 * d$x) + 1.5 * (1 - exp(-0.6 * d$z)) + 0.5 * exp(0.3 * d$w)
  d$y <- switch(fam, gaussian = lat + rnorm(n, 0, 0.5),
                poisson = rpois(n, pmax(lat, 0.05)),
                bernoulli = rbinom(n, 1, plogis(lat - 2)),
                cumulative = ordered(cut(lat + rlogis(n),
                                         c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4)))
  d
}
cat("ROWHEAD body fam seed in_order reversed renamed renamed_reversed maxdiff\n")
for (bn in c("six", "five", "tie")) {
  terms <- switch(bn, six = terms6, five = terms5, tie = termsT)
  np <- if (bn == "five") 5 else 6
  for (fam in names(F)) for (seed in c(501, 502, 503)) {
    d <- gen(fam, seed, tie = bn == "tie")
    v <- c()
    for (sp in list(list("b", FALSE), list("b", TRUE), list("rename", FALSE),
                    list("rename", TRUE))) {
      s <- spell(terms, np, maps[[sp[[1]]]], sp[[2]])
      fm <- bf(as.formula(paste("y ~", s$body)),
               as.formula(paste(s$decl, "~ 1")), nl = TRUE)
      v <- c(v, llv(tt(frm(fm, family = F[[fam]], data = d))))
    }
    cat(sprintf("ROW %s %s %d %.6f %.6f %.6f %.6f %.3g\n", bn, fam, seed,
                v[1], v[2], v[3], v[4], max(v) - min(v)))
    flush.console()
  }
}
cat("DONE\n")
