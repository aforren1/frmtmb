# Lane arspat: what brms 2.23.0 writes for each combination the lane
# looks at. make_stancode() only; no compile. Output:
# dev/arspat-out/brms-probe.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(2)
g <- rep(1:6, each = 8)
t <- unlist(lapply(1:6, function(i) if (i %% 2) c(1:4, 6:9) else 1:8))
d <- data.frame(g = factor(g), t = t, x = rnorm(48))
d$y <- d$x + rnorm(48)
d$s <- runif(48, 0.2, 0.5)
d$ym <- d$y; d$ym[c(3, 20)] <- NA
show <- function(label, expr, pat = NULL) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(expr, error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n"); NULL
  })
  if (is.null(r)) return(invisible(NULL))
  code <- as.character(r)
  lines <- strsplit(code, "\n")[[1]]
  keep <- if (is.null(pat)) grep("lpdf|target|sigma|se2|Yl|Ymi|cor|ar\\b",
                                  lines, value = TRUE) else
    grep(pat, lines, value = TRUE)
  cat(utils::tail(keep, 30), sep = "\n")
  invisible(r)
}
show("se(sigma=TRUE) + ar cov=TRUE",
     make_stancode(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE),
                   data = d))
show("se(sigma=FALSE) + ar cov=TRUE",
     make_stancode(y | se(s) ~ x + ar(t, g, cov = TRUE), data = d))
show("se + ma cov=TRUE",
     make_stancode(y | se(s, sigma = TRUE) ~ x + ma(t, g, cov = TRUE),
                   data = d))
show("se + cosy",
     make_stancode(y | se(s, sigma = TRUE) ~ x + cosy(t, g), data = d))
show("se + unstr",
     make_stancode(y | se(s, sigma = TRUE) ~ x + unstr(t, g), data = d))
show("se + ar cov=TRUE, sigma ~ x",
     make_stancode(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE),
                      sigma ~ x), data = d))
show("student se + ar cov=TRUE",
     make_stancode(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE),
                   data = d, family = student()))
show("mi + ar cov=TRUE",
     make_stancode(ym | mi() ~ x + ar(t, g, cov = TRUE), data = d))
show("mi + se + ar cov=TRUE",
     make_stancode(ym | mi(s) ~ x + ar(t, g, cov = TRUE), data = d))
show("mi + cosy",
     make_stancode(ym | mi() ~ x + cosy(t, g), data = d))
show("student mi + ar cov=TRUE",
     make_stancode(ym | mi() ~ x + ar(t, g, cov = TRUE), data = d,
                   family = student()))
# sar
W <- matrix(0, 48, 48)
for (i in 1:47) { W[i, i + 1] <- 1; W[i + 1, i] <- 1 }
W <- W / rowSums(W)
show("sar error + mi",
     make_stancode(ym | mi() ~ x + sar(W, type = "error"), data = d,
                   data2 = list(W = W)))
show("sar lag + mi",
     make_stancode(ym | mi() ~ x + sar(W, type = "lag"), data = d,
                   data2 = list(W = W)))
show("sar error + ar cov=FALSE",
     make_stancode(y ~ x + ar(t, g) + sar(W, type = "error"), data = d,
                   data2 = list(W = W)))
show("sar lag + ar cov=FALSE",
     make_stancode(y ~ x + ar(t, g) + sar(W, type = "lag"), data = d,
                   data2 = list(W = W)))
show("sar + ar cov=TRUE",
     make_stancode(y ~ x + ar(t, g, cov = TRUE) + sar(W), data = d,
                   data2 = list(W = W)))
show("sar error + se",
     make_stancode(y | se(s, sigma = TRUE) ~ x + sar(W, type = "error"),
                   data = d, data2 = list(W = W)))
show("fcor + se",
     make_stancode(y | se(s, sigma = TRUE) ~ x + fcor(M), data = d,
                   data2 = list(M = diag(48))))
show("fcor + mi",
     make_stancode(ym | mi() ~ x + fcor(M), data = d,
                   data2 = list(M = diag(48))))
