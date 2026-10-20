# Lane prilabel: which per-column-sd structures take (1 + x | g).
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(sample(1:20, 200, TRUE)))
dd$y <- rnorm(200)
for (cs in c("cs", "toep", "hetar1", "homcs", "ar1")) {
  f <- as.formula(paste0("y ~ x + ", cs, "(1 + x | g)"))
  r <- tryCatch({
    des <- frmtmb:::prior_design(bf(f), dd, gaussian(), list())
    bk <- des$frame$re_blocks[[1]]
    paste(bk$covstruct, "sd coefs:",
          paste(frmtmb:::block_sd_coefs(des$frame, bk), collapse = ","))
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat(cs, ":", r, "\n")
}
