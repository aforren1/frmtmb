# Spot check, part 3: the prior interface, which the plan puts in bin 1.
# brms 2.23 spells it default_prior(); get_prior() is the alias kept
# since 2.20.14. Everything below is brms's own assertion with the data
# swapped for one frmtmb can build a frame from.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
load("dev/brms-suite/brms/data/epilepsy.rda")

say <- function(label, expr) {
  out <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("---", label, "---\n")
  print(out)
  cat("\n")
}

cat("## names brms exports that frmtmb does not\n")
for (nm in c("default_prior", "validate_prior", "get_prior", "set_prior",
             "prior", "prior_", "prior_string", "prior_summary",
             "empty_prior", "as.brmsprior")) {
  cat(sprintf("  %-16s in frmtmb: %s\n", nm,
              nm %in% getNamespaceExports("frmtmb")))
}

cat("\n## brms: default_prior() finds all classes for which priors can be set\n")
say("get_prior(count ~ zBase * Trt + (1|patient) + (1+Trt|visit), poisson)", {
  p <- get_prior(count ~ zBase * Trt + (1 | patient) + (1 + Trt | visit),
                 data = epilepsy, family = poisson())
  list(class_of_result = class(p), columns = names(p),
       classes = sort(unique(p$class)))
})
cat("brms's answer for the same call, from tests.priors.R:\n")
cat("  sort(c(rep('b', 4), 'cor', 'cor', 'Intercept', rep('sd', 6)))\n\n")

cat("## brms: set_prior allows arguments to be vectors\n")
say("set_prior('normal(0, 2)', class = c('b', 'sd'))", {
  bp <- set_prior("normal(0, 2)", class = c("b", "sd"))
  list(class = class(bp), prior = bp$prior, cls = bp$class)
})

cat("## brms: print for class brmsprior\n")
say("print(set_prior('normal(0,1)'))", capture.output(print(set_prior("normal(0,1)"))))
say("print(set_prior('normal(0,1)', coef = 'x'))",
    capture.output(print(set_prior("normal(0,1)", coef = "x"))))

cat("## brms: set_prior alias functions produce equivalent results\n")
say("set_prior(class='sd') == prior(normal(0, 1), class = sd)",
    all.equal(set_prior("normal(0, 1)", class = "sd"),
              prior(normal(0, 1), class = sd)))
say("set_prior(...) == prior_string(...)",
    all.equal(set_prior("normal(0, 1)", class = "sd"),
              prior_string("normal(0, 1)", class = "sd")))
say("set_prior(nlpar) == prior_(~normal(0, 1), class = ~sd, nlpar = quote(a))",
    tryCatch(all.equal(set_prior("normal(0, 1)", class = "sd", nlpar = "a"),
                       prior_(~normal(0, 1), class = ~sd, nlpar = quote(a))),
             error = function(e) conditionMessage(e)))

cat("## brms: prior + prior concatenates\n")
say("prior(normal(0,10), class = b) + prior(cauchy(0,2), class = sd)",
    tryCatch({
      p <- prior(normal(0, 10), class = b) + prior(cauchy(0, 2), class = sd)
      list(class = class(p), n = NROW(p))
    }, error = function(e) conditionMessage(e)))
