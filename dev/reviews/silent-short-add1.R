# Short check, MAJOR 2 and MINOR 3: add1() bar spellings, dharma args.
.libPaths(c("C:/Users/adf44/source/r/wt-silent-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
ll <- function(f) as.numeric(logLik(f))
dfl <- function(f) attr(logLik(f), "df")
set.seed(11)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)),
                g2 = factor(sample(1:3, n, TRUE)),
                f = factor(sample(c("a", "b"), n, TRUE)))
re <- rnorm(8)[d$g] + rnorm(24, 0, 0.5)[interaction(d$g, d$g2)]
d$y <- 1 + 0.5 * d$x + re + 0.4 * rnorm(8)[d$g] * d$x + rnorm(n, 0, 0.5)
f1 <- frm(y ~ x, data = d)
cases <- list(
  "(x | g)" = y ~ x + (x | g),
  "(1 | g/g2)" = y ~ x + (1 | g/g2),
  "(1 | g:g2)" = y ~ x + (1 | g:g2),
  "(1 | gr(g, by = f))" = y ~ x + (1 | gr(g, by = f)),
  "(x || g)" = y ~ x + (x || g)
)
for (nm in names(cases)) {
  sc <- stats::as.formula(paste("~ . +", nm))
  h <- tryCatch(frm(cases[[nm]], data = d), error = function(e) e)
  if (inherits(h, "error")) {
    cat(sprintf("%-22s hand refit ERROR %s\n", nm,
                substr(conditionMessage(h), 1, 100)))
    next
  }
  for (form in c("formula", "character")) {
    s <- if (form == "formula") sc else {
      # step()'s spelling: the bare label add.scope() makes
      attr(stats::terms(stats::update.formula(y ~ x, sc)), "term.labels")[-1]
    }
    a <- tryCatch(suppressWarnings(add1(f1, s, test = "Chisq")),
                  error = function(e) e)
    if (inherits(a, "error")) {
      cat(sprintf("%-22s %-9s ERROR %s\n", nm, form,
                  substr(conditionMessage(a), 1, 120)))
      next
    }
    cat(sprintf(paste0("%-22s %-9s rows [%s] Df %s LRT %.4f | hand Df %d",
                       " LRT %.4f\n"), nm, form,
                paste(rownames(a)[-1], collapse = "; "),
                paste(a$Df[-1], collapse = ","), a$LRT[2],
                dfl(h) - dfl(f1), 2 * (ll(h) - ll(f1))))
  }
}
# a fixed term holding a logical or
d$b1 <- d$x > 0; d$b2 <- d$z > 0
a <- tryCatch(add1(f1, ~ . + I(b1 | b2), test = "Chisq"),
              error = function(e) conditionMessage(e))
h <- frm(y ~ x + I(b1 | b2), data = d)
if (is.character(a)) cat("I(b1 | b2) ERROR", a, "\n") else
  cat(sprintf("I(b1 | b2) rows [%s] LRT %.4f hand %.4f\n",
              paste(rownames(a)[-1], collapse = "; "), a$LRT[2],
              2 * (ll(h) - ll(f1))))
# step() forward with a bar term in its upper scope
st <- tryCatch(step(f1, scope = list(lower = ~ x, upper = ~ x + z + (1 | g)),
                    direction = "forward", trace = 0),
               error = function(e) conditionMessage(e))
cat("STEP forward:", if (is.character(st)) st else
  deparse(st$bform$formula), "\n")
# drop1 with a bar term as scope
cat("DROP1 bar scope:", tryCatch({drop1(frm(y ~ x + (1 | g), data = d),
                                        scope = "1 | g"); "ran"},
                                 error = function(e) conditionMessage(e)), "\n")

## MINOR 3: dharma args
set.seed(5)
xx <- rnorm(300); ys <- 1 + 0.5 * xx + rnorm(300)
ct <- runif(300, 0.8, 2.5)
cen <- as.numeric(ys > ct)
dc <- data.frame(x = xx, y = ifelse(cen == 1, ct, ys), cen = cen)
fc <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = dc)
dn <- dc; dn$cen <- 0; dn$y <- ys
fn <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = dn)
fp <- frm(y ~ x, data = dn)
try_d <- function(f, ...) tryCatch({
  r <- suppressMessages(dharma_residuals(f, nsim = 100, seed = 1, ...))
  "ran"
}, error = function(e) paste("refused:", substr(conditionMessage(e), 1, 70)))
for (a in list(list(method = "traditional"), list(method = "trad"),
               list(method = "t"), list(method = "PIT"),
               list(method = c("PIT", "traditional")),
               list(rotation = "estimated"))) {
  lab <- paste(names(a), unlist(a), collapse = " ")
  cat(sprintf("DHARMA %-32s censored: %-40s | cens() none censored: %s | no cens(): %s\n",
              lab, do.call(try_d, c(list(fc), a)),
              do.call(try_d, c(list(fn), a)), do.call(try_d, c(list(fp), a))))
}
