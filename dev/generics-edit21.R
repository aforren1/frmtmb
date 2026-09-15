root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
sub1 <- function(path, old, new) {
  f <- file.path(root, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (length(gregexpr(old, txt, fixed = TRUE)[[1]]) != 1L)
    stop("many in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("edited", path, "\n")
}

# tickof() reported NA for proc.time(): 4000 reads of a 20 ms clock can
# finish inside one tick, so there were no positive differences to
# take a minimum of. It now samples until it has seen the clock move
# enough times, which is the point of measuring a tick at all.
sub1("dev/generics-loadcost.R",
paste0("tickof <- function(f) {\n",
       "  v <- replicate(4000, f())\n",
       "  d <- diff(v)\n",
       "  d <- d[d > 0]\n",
       "  if (!length(d)) return(c(min = NA_real_, median = NA_real_))\n",
       "  c(min = min(d), median = stats::median(d))\n",
       "}\n"),
paste0("tickof <- function(f, want = 200L, budget = 5) {\n",
       "  d <- numeric(0)\n",
       "  t0 <- Sys.time()\n",
       "  repeat {\n",
       "    v <- replicate(20000, f())\n",
       "    e <- diff(v)\n",
       "    d <- c(d, e[e > 0])\n",
       "    if (length(d) >= want) break\n",
       "    if (as.numeric(difftime(Sys.time(), t0,\n",
       "                            units = \"secs\")) > budget) break\n",
       "  }\n",
       "  if (!length(d)) return(c(min = NA_real_, median = NA_real_))\n",
       "  c(min = min(d), median = stats::median(d))\n",
       "}\n"))
cat("DONE\n")
