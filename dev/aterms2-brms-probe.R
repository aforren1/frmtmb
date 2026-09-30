.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
try_show <- function(label, expr) {
  cat("\n==========", label, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  print(r)
  invisible(r)
}
set.seed(1)
n <- 12
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), x = rnorm(n), z = rnorm(n),
                g = factor(rep(letters[1:4], 3)),
                f = factor(rep(c("u", "v", "w"), each = 4)),
                s1 = rep(c(TRUE, FALSE), 6), s2 = c(rep(TRUE, 8), rep(FALSE, 4)))
d$y2[!d$s2] <- NA
d$z[!d$s2] <- NA

# --- subset, basic -----------------------------------------------------
bf1 <- bf(y1 | subset(s1) ~ x + (1 | g)) + bf(y2 | subset(s2) ~ z + (1 | g)) +
  set_rescor(FALSE)
sd1 <- try_show("mv subset standata", standata(bf1, d))
try_show("names", names(sd1))
try_show("N_1 J_1_y1 J_1_y2 N_y1 N_y2",
         sd1[c("N", "N_1", "J_1_y1", "J_2_y2", "N_2", "N_y1", "N_y2")])
cat(stancode(bf1, d))

# a level of g present only in rows excluded by s2 (g = d rows 4, 8, 12?)
d2 <- d
d2$g <- factor(ifelse(!d2$s2, "only_out", as.character(d2$g)))
bf2 <- bf(y1 ~ x) + bf(y2 | subset(s2) ~ z + (1 | g)) + set_rescor(FALSE)
sd2 <- try_show("level only outside subset", standata(bf2, d2))
try_show("N_1 levels, J_1_y2", sd2[c("N_1", "J_1_y2")])

# the same grouping term in two subsetted responses with |ID|
bf2b <- bf(y1 | subset(s1) ~ x + (1 | p | g)) +
  bf(y2 | subset(s2) ~ z + (1 | p | g)) + set_rescor(FALSE)
sd2b <- try_show("|ID| with subsets", standata(bf2b, d2))
try_show("N_1 J_1_y1 J_1_y2", sd2b[c("N_1", "J_1_y1", "J_1_y2")])

# fixed-effect factor with an empty level in a subset
bf3 <- bf(y1 ~ x) + bf(y2 | subset(s2) ~ f) + set_rescor(FALSE)
sd3 <- try_show("factor level empty in subset", standata(bf3, d))
try_show("X_y2", sd3$X_y2)

# univariate subset
try_show("univariate subset", standata(bf(y1 | subset(s1) ~ x), d))
try_show("univariate subset N", standata(bf(y1 | subset(s1) ~ x), d)$N)

# rescor with subset
try_show("rescor subset", standata(bf(y1 | subset(s1) ~ x) +
                                     bf(y2 | subset(s2) ~ z) +
                                     set_rescor(TRUE), d))
# NA in y1 outside its own subset
d4 <- d; d4$y1[2] <- NA
try_show("NA in y1 on an excluded row",
         standata(bf1, d4)$N_y1)
# non-logical subset
d5 <- d; d5$s1 <- as.numeric(d5$s1)
try_show("numeric subset", standata(bf1, d5)$N_y1)
d5$s1[1] <- 2
try_show("numeric subset with 2", standata(bf1, d5)$N_y1)
# subset with an expression
try_show("subset expression", standata(bf(y1 | subset(x > 0) ~ 1) +
                                        bf(y2 | subset(s2) ~ z) +
                                        set_rescor(FALSE), d)$N_y1)

# --- mi idx -----------------------------------------------------------
set.seed(2)
dm <- data.frame(y = rnorm(10), x = c(rnorm(9), NA), z = rnorm(10),
                 g1 = sample(1:5, 10, TRUE), g2 = 10:1, g3 = 1:10,
                 s = c(FALSE, rep(TRUE, 9)))
bm <- bf(y ~ mi(x, idx = g1)) + bf(x | mi() + index(g2) + subset(s) ~ 1) +
  set_rescor(FALSE)
sdm <- try_show("mi idx standata", standata(bm, dm))
try_show("mi idx parts", sdm[grep("idx|Jmi|Nmi|N_|Y_", names(sdm))])
cat(stancode(bm, dm))
# index without subset
bm2 <- bf(y ~ mi(x, idx = g1)) + bf(x | mi() + index(g3) ~ 1) +
  set_rescor(FALSE)
sdm2 <- try_show("mi idx no subset", standata(bm2, dm))
try_show("idxl", sdm2[grep("idx", names(sdm2))])
# idx values that do not match
dm3 <- dm; dm3$g1[1] <- 99
try_show("unmatched idx", standata(bm2, dm3))
# duplicated index values
dm4 <- dm; dm4$g3[2] <- 1
try_show("duplicated index", standata(bm2, dm4))
# index on a response that no mi() uses
try_show("index unused", standata(bf(y ~ x) + bf(x | index(g3) ~ 1) +
                                   set_rescor(FALSE), dm[-10, ])$N)
# idx on a univariate mi
try_show("mi(x, idx) univariate, x not response",
         standata(bf(y ~ mi(x, idx = g1)), dm))
# subset on the mi() user with idx
bm5 <- bf(y | subset(s) ~ mi(x, idx = g1)) + bf(x | mi() + index(g3) ~ 1) +
  set_rescor(FALSE)
sdm5 <- try_show("user subsetted with idx", standata(bm5, dm))
try_show("idxl", sdm5[grep("idx|N_", names(sdm5))])
# mi(x, idx) with a response variable x having mi and no NA
# me() in a subsetted formula
try_show("me in subsetted",
         standata(bf(y | subset(s) ~ me(z, 1)) + bf(x ~ 1) +
                    set_rescor(FALSE), dm[-10, ]))

# --- rate ---------------------------------------------------------------
dr <- data.frame(y = rpois(10, 1), x = rnorm(10), time = 1:10,
                 cc = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0))
cat("\n=== rate poisson\n")
cat(stancode(y | rate(time) ~ x, dr, poisson()))
cat("\n=== rate poisson identity\n")
cat(stancode(y | rate(time) ~ x, dr, poisson("identity")))
cat("\n=== rate negbinomial\n")
cat(stancode(y | rate(time) ~ x, dr, negbinomial()))
cat("\n=== rate negbinomial shape ~ x\n")
cat(stancode(bf(y | rate(time) ~ x, shape ~ x), dr, negbinomial()))
cat("\n=== rate negbinomial2\n")
cat(stancode(y | rate(time) ~ x, dr, negbinomial2()))
cat("\n=== rate geometric\n")
cat(stancode(y | rate(time) ~ x, dr, geometric()))
cat("\n=== rate poisson cens\n")
cat(stancode(y | rate(time) + cens(cc) ~ x, dr, poisson()))
try_show("rate gaussian", stancode(y | rate(time) ~ x, dr, gaussian()))
try_show("rate zip", stancode(y | rate(time) ~ x, dr, zero_inflated_poisson()))
try_show("rate 0", standata(y | rate(time - 1) ~ x, dr, poisson()))
try_show("rate const", standata(y | rate(2) ~ x, dr, poisson())$denom)

# --- cat --------------------------------------------------------------
dc <- data.frame(s = sample(1:5, 9, TRUE), x = rnorm(9))
try_show("cat standata", {
  w <- NULL
  r <- withCallingHandlers(standata(s | cat(6) ~ x, dc, cumulative()),
                           warning = function(ww) {
                             w <<- c(w, conditionMessage(ww))
                             invokeRestart("muffleWarning")
                           })
  list(nthres = r$nthres, warnings = w)
})
try_show("cat on gaussian", standata(s | cat(6) ~ x, dc, gaussian()))
try_show("cat + thres", standata(s | cat(6) + thres(5) ~ x, dc, cumulative()))
