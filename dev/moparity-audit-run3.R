# The per-scale floor's audit, narrowed (dev/moparity-findings.md, item
# 5): as dev/moparity-audit-run.R, but recording the smallest theta at
# the positions the floor changes, the log standard deviations of the
# blocks that read sd_floored() (us and us_t with two or more
# coefficients, cs, toep, ar1, hetar1, gr_cov and gr_prec with two or
# more). The blocks are found in the environment of the function
# MakeADFun() tapes, where build_objective() keeps them.
#   Rscript dev/moparity-audit-run2.R <package> <test file> <out file>
local({
LIB <- "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER"))) Sys.setenv(R_MAKEVARS_USER = mk)
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
out <- a[3]
aud <- new.env()
aud$min <- Inf
aud$min_all <- Inf
aud$n <- 0
aud$n_found <- 0
aud$n_theta_obj <- 0
aud$hits <- character(0)
find_blocks <- function(fn, depth = 0L, seen = list()) {
  if (!is.function(fn) || depth > 4L) return(NULL)
  e <- environment(fn)
  if (is.null(e) || isNamespace(e)) return(NULL)
  if (exists("blocks", envir = e, inherits = FALSE) &&
        exists("block_fns", envir = e, inherits = FALSE)) {
    return(get("blocks", envir = e))
  }
  for (nm in ls(e, all.names = TRUE)) {
    v <- tryCatch(get(nm, envir = e), error = function(err) NULL)
    if (is.function(v) && !identical(v, fn)) {
      r <- find_blocks(v, depth + 1L)
      if (!is.null(r)) return(r)
    }
  }
  NULL
}
floored_pos <- function(blocks) {
  reg <- get("covstruct_registry", asNamespace("frmtmb"))
  out <- integer(0)
  for (bk in blocks) {
    cs <- bk[["covstruct"]]
    d <- bk[["dim"]]
    hit <- (cs %in% c("us", "us_t", "gr_cov", "gr_prec") && d >= 2L) ||
      cs %in% c("cs", "toep", "ar1", "hetar1")
    if (hit) out <- c(out, bk[["theta_idx"]][reg[[cs]]$sd_idx(d)])
  }
  out
}
orig <- RTMB::MakeADFun
wrapper <- function(...) {
  obj <- orig(...)
  pn <- names(obj$par)
  th <- which(pn == "theta")
  if (!length(th)) return(obj)
  aud$n_theta_obj <- aud$n_theta_obj + 1
  args <- list(...)
  blocks <- tryCatch(find_blocks(args[[1]]), error = function(e) NULL)
  fp <- if (length(blocks)) th[floored_pos(blocks)] else integer(0)
  fp <- fp[!is.na(fp)]
  if (length(blocks)) aud$n_found <- aud$n_found + 1
  for (nm in c("fn", "gr", "he")) {
    g <- obj[[nm]]
    if (!is.function(g)) next
    obj[[nm]] <- local({
      g0 <- g
      function(...) {
        args <- list(...)
        x <- if (length(args)) args[[1]] else NULL
        if (is.numeric(x) && length(x) == length(pn)) {
          aud$min_all <- min(aud$min_all, x[th])
          if (length(fp)) {
            aud$min <- min(aud$min, x[fp])
            if (min(x[fp]) < -308.45) {
              aud$hits <- unique(c(aud$hits, paste(vapply(blocks, `[[`, "", "term_label"), collapse = " + ")))
            }
          }
          aud$n <- aud$n + 1
        }
        g0(...)
      }
    })
  }
  obj
}
utils::assignInNamespace("MakeADFun", wrapper, "RTMB")
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
res <- tryCatch(as.data.frame(test_file(f, package = p,
                                        env = testthat::test_env(p),
                                        reporter = "silent")),
                error = function(e) NULL)
cat(sprintf(paste("AUDIT2 %s %s calls %d objects %d found %d",
                  "min_floored %s min_theta %s fail %s\n"),
            p, basename(f), aud$n, aud$n_theta_obj, aud$n_found,
            format(aud$min, digits = 10), format(aud$min_all, digits = 10),
            if (is.null(res)) "LOADERROR" else sum(res$failed)),
    file = out)
cat(paste0("HIT ", aud$hits, "
"), file = out, append = TRUE)
})
