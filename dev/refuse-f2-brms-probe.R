# What brms 2.23.0 says on the two id questions of fix round 2. Seed 11800.
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
set.seed(11800)
n <- 120
d <- data.frame(y = rnorm(n), x = rnorm(n), w1 = runif(n), w2 = runif(n),
                g = sample(letters[1:6], n, TRUE),
                h = sample(LETTERS[1:6], n, TRUE),
                g1 = sample(letters[1:6], n, TRUE),
                g2 = sample(letters[1:6], n, TRUE))
d$W <- cbind(d$w1, d$w2)
try_code <- function(lab, f) {
  r <- tryCatch({
    code <- brms::make_stancode(f, data = d)
    sprintf("builds; cor_1 %s", grepl("cor_1", code))
  }, error = function(e) paste("refuses:", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", lab, substr(r, 1, 200)))
}
try_code("cor FALSE over another", brms::bf(
  y ~ x + (1 | gr(g, id = "p")), sigma ~ (1 | gr(h, id = "p", cor = FALSE))))
try_code("cor FALSE first", brms::bf(
  y ~ x + (1 | gr(g, id = "p", cor = FALSE)), sigma ~ (1 | gr(h, id = "p"))))
try_code("both cor FALSE", brms::bf(
  y ~ x + (1 | gr(g, id = "p", cor = FALSE)),
  sigma ~ (1 | gr(h, id = "p", cor = FALSE))))
try_code("|p| spelling", brms::bf(
  y ~ x + (1 | p | gr(g, cor = FALSE)), sigma ~ (1 | p | gr(h, cor = FALSE))))
try_code("W beside cbind", brms::bf(
  y ~ x + (1 | mm(g1, g2, id = "p", weights = W)),
  sigma ~ (1 | mm(g1, g2, id = "p", weights = cbind(w1, w2)))))
try_code("different weights", brms::bf(
  y ~ x + (1 | mm(g1, g2, id = "p", weights = cbind(w1, w2))),
  sigma ~ (1 | mm(g1, g2, id = "p", weights = cbind(w2, w1)))))
