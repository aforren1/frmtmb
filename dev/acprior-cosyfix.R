# Fix round 1, item 2: frm_sample()'s cosy near its bound at 0, after
# the bound chart, against brms 2.23.0 on the same data and seeds. The
# reviewer's construction (dev/review-acprior-cosydiv.R): the cosy data
# of dev/acprior-agree.R, 4 chains of 2000, seeds 3 to 5, mode and
# random starts. Arg 1: "frm" (lane library) or "brms".
# Log: dev/acprior-out/cosyfix-<arm>.txt
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
src <- readLines(file.path(wt, "dev", "acprior-agree.R"))
i0 <- grep("^make_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
dc <- make_data()$cosy
fo <- y ~ x + cosy(time = t, gr = g)
if (identical(arm, "frm")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
  fit <- frm(fo, data = dc)
  for (s in 3:5) {
    for (init in c("last.par.best", "random")) {
      w <- character(0)
      ds <- withCallingHandlers(
        suppressMessages(frm_sample(fit, chains = 4, iter = 2000, seed = s,
                                    refresh = 0, init = init)),
        warning = function(x) {
          w <<- c(w, conditionMessage(x))
          invokeRestart("muffleWarning")
        })
      nd <- sum(sapply(rstan::get_sampler_params(ds$stanfit,
                                                 inc_warmup = FALSE),
                       function(m) sum(m[, "divergent__"])))
      cz <- as.matrix(ds)[, "cosy"]
      cat(sprintf("init %-13s seed %d: %d divergent; cosy mean %.4f sd %.4f min %.2e\n",
                  init, s, nd, mean(cz), sd(cz), min(cz)))
      if (length(w)) cat("  warnings:", paste(w, collapse = " || "), "\n")
    }
  }
} else {
  suppressPackageStartupMessages(library(brms))
  code <- stancode(fo, data = dc)
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  path <- file.path(Sys.getenv("FRMTMB_STAN_CACHE"),
                    paste0(unname(tools::md5sum(f)), ".rds"))
  mod <- if (file.exists(path)) readRDS(path) else {
    mm <- rstan::stan_model(model_code = code, save_dso = TRUE)
    saveRDS(mm, path)
    mm
  }
  for (s in 3:5) {
    sf <- suppressWarnings(rstan::sampling(mod, data = standata(fo, data = dc),
                                           chains = 4, iter = 2000, seed = s,
                                           refresh = 0))
    cz <- as.matrix(sf, pars = "cosy")[, 1]
    cat(sprintf("brms seed %d: %d divergent; cosy mean %.4f sd %.4f\n", s,
                rstan::get_num_divergent(sf), mean(cz), sd(cz)))
  }
}
