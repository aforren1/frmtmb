# Reviewer round 2: do the lane's density tests catch a small relative
# error? Runs the three density blocks of test-simplex-families.R with
# each family's lpdf multiplied by (1 + eps), recording every
# expect_lt() value. Output: dev/famsimplex-rev2-ulps.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
tf <- readLines(paste0("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/",
                       "tests/testthat/test-simplex-families.R"))
start <- grep("^test_that\\(\"the dirichlet density", tf)
end <- grep("^test_that\\(\"a dirichlet fit's logLik", tf) - 1L
helpers <- tf[seq_len(grep("^sx_data_dir", tf) - 1L)]
blocks <- tf[start:end]
for (f in c("fam_dirichlet_impl", "fam_dirichlet2_impl",
            "fam_logistic_normal_impl")) {
  blocks <- gsub(paste0("frmtmb:::", f, "("),
                 paste0("pert(frmtmb:::", f, ")("), blocks, fixed = TRUE)
}
cat("blocks:", length(grep("^test_that", blocks)), " pert calls:",
    length(grep("pert(", blocks, fixed = TRUE)), "\n")
for (eps in c(0, 1e-13, 1e-12, 1e-11)) {
  env <- new.env(parent = asNamespace("frmtmb"))
  rec <- list()
  env$pert <- function(ctor) {
    function(...) {
      fam <- ctor(...)
      lp <- fam$lpdf
      fam$lpdf <- function(...) lp(...) * (1 + eps)
      fam
    }
  }
  env$test_that <- function(desc, code) {
    cur <<- desc
    invisible(code)
  }
  cur <- ""
  env$expect_lt <- function(object, expected, label = NULL, ...) {
    rec[[length(rec) + 1L]] <<- data.frame(test = substr(cur, 5, 22),
                                           value = object, bound = expected)
  }
  eval(parse(text = helpers), env)
  eval(parse(text = blocks), env)
  r <- do.call(rbind, rec)
  cat(sprintf("\neps %g: %d expectations, %d over the bound\n", eps,
              nrow(r), sum(r$value >= r$bound)))
  print(stats::aggregate(value ~ test, r, max))
}
