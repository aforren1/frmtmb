# Lane meme: what brms 2.23.0 does with me() terms. No compile:
# default_prior(), stancode() and standata() only.
# Output: dev/meme-out/brms.txt (Rscript dev/meme-brms.R > that file).
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(11)
n <- 40
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                sx = 0.3, sz = runif(n, 0.2, 0.4),
                g = factor(rep(1:8, each = 5)), h = factor(rep(1:5, 8)))
# grouped me() needs x and sdx constant within a level
d$xg <- rnorm(8)[d$g]
d$zg <- rnorm(8)[d$g]
d$sg <- 0.25
d$y <- rnorm(n)

show <- function(lab, f, prior = NULL, sample_prior = "no") {
  cat("\n==================== ", lab, "\n", sep = "")
  dp <- tryCatch(default_prior(f, data = d), error = function(e) e)
  if (inherits(dp, "error")) {
    cat("default_prior ERROR:", conditionMessage(dp), "\n")
    return(invisible())
  }
  print(as.data.frame(dp)[, c("prior", "class", "coef", "group", "lb",
                              "ub", "source")])
  sc <- tryCatch(stancode(f, data = d, prior = prior,
                          sample_prior = sample_prior),
                 error = function(e) e)
  if (inherits(sc, "error")) {
    cat("stancode ERROR:", conditionMessage(sc), "\n")
    return(invisible())
  }
  lines <- strsplit(as.character(sc), "\n")[[1]]
  keep <- grep("me|Xn|noise|Lme|prior_", lines, value = TRUE)
  cat("-- stancode lines mentioning me/Xn/noise/Lme/prior_:\n")
  cat(keep, sep = "\n")
  sd <- standata(f, data = d, prior = prior, sample_prior = sample_prior)
  nm <- grep("me|Xn|noise", names(sd), value = TRUE)
  cat("-- standata me names:", nm, "\n")
  for (k in nm) cat("   ", k, ":", head(sd[[k]], 6), "\n")
  invisible(sc)
}

show("A one me()", bf(y ~ me(x, sx)))
show("B two me(), correlated", bf(y ~ me(x, sx) + me(z, sz)))
show("C two me(), set_mecor(FALSE)",
     bf(y ~ me(x, sx) + me(z, sz)) + set_mecor(FALSE))
show("D one me(gr = g)", bf(y ~ me(xg, sg, gr = g)))
show("E two me(gr = g), correlated",
     bf(y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g)))
show("F me() + (1 | h)", bf(y ~ me(x, sx) + (1 | h)))
show("G me() and me(gr = g): two groups",
     bf(y ~ me(x, sx) + me(z, sz) + me(xg, sg, gr = g)))
show("H interaction me(x):me(z)", bf(y ~ me(x, sx) * me(z, sz)))

cat("\n==================== I missing sdx\n")
r <- tryCatch(default_prior(bf(y ~ me(x)), data = d),
              error = function(e) conditionMessage(e))
print(r)
d2 <- d
d2$sx <- c(NA, rep(0.3, n - 1))
cat("-- NA in sdx:\n")
r <- tryCatch({
  s <- standata(bf(y ~ me(x, sx)), data = d2)
  c(N = s$N, noise1 = s$noise_1[1:3])
}, error = function(e) conditionMessage(e))
print(r)
d3 <- d
d3$x[1] <- NA
cat("-- NA in x:\n")
r <- tryCatch({
  s <- standata(bf(y ~ me(x, sx)), data = d3)
  c(N = s$N, Xn = s$Xn_1[1:3])
}, error = function(e) conditionMessage(e))
print(r)
cat("-- sdx = 0 in a row:\n")
d4 <- d
d4$sx[2] <- 0
r <- tryCatch({
  s <- standata(bf(y ~ me(x, sx)), data = d4)
  c(N = s$N, noise = s$noise_1[1:3])
}, error = function(e) conditionMessage(e))
print(r)

cat("\n==================== J user priors, per coefficient\n")
pr <- c(set_prior("normal(0, 2)", class = "meanme"),
        set_prior("normal(1, 0.5)", class = "meanme", coef = "mez"),
        set_prior("exponential(1)", class = "sdme"),
        set_prior("student_t(3, 0, 1)", class = "sdme", coef = "mex"),
        set_prior("lkj(2)", class = "corme"))
show("J priors on B", bf(y ~ me(x, sx) + me(z, sz)), prior = pr)
cat("-- validate_prior:\n")
print(as.data.frame(validate_prior(pr, bf(y ~ me(x, sx) + me(z, sz)),
                                   data = d))[, c("prior", "class", "coef",
                                                  "group", "lb", "ub",
                                                  "source")])
cat("\n-- corme with group on E:\n")
pr2 <- c(set_prior("lkj(3)", class = "corme", group = "g"),
         set_prior("normal(0, 1)", class = "meanme", coef = "mexg"))
show("K priors on E", bf(y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g)),
     prior = pr2)
cat("\n-- a bound on sdme:\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx)), data = d,
                       prior = set_prior("normal(0, 1)", class = "sdme",
                                         ub = 3)),
              error = function(e) conditionMessage(e))
cat(grep("sdme", strsplit(as.character(r), "\n")[[1]], value = TRUE),
    sep = "\n")
cat("\n-- lb on meanme:\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx)), data = d,
                       prior = set_prior("normal(0, 1)", class = "meanme",
                                         lb = 0)),
              error = function(e) conditionMessage(e))
cat(grep("meanme", strsplit(as.character(r), "\n")[[1]], value = TRUE),
    sep = "\n")
cat("\n-- coef-level bound on sdme:\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                       prior = set_prior("normal(0, 1)", class = "sdme",
                                         coef = "mex", ub = 3)),
              error = function(e) conditionMessage(e))
print(if (is.character(r) && length(r) == 1L && !inherits(r, "brmsmodel"))
  r else grep("sdme", strsplit(as.character(r), "\n")[[1]], value = TRUE))
cat("\n-- corme coef:\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                       prior = set_prior("lkj(2)", class = "corme",
                                         coef = "mex")),
              error = function(e) conditionMessage(e))
print(r)
cat("\n-- non-lkj on corme:\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                       prior = set_prior("normal(0, 1)", class = "corme")),
              error = function(e) conditionMessage(e))
print(substr(paste(r, collapse = " "), 1, 300))
cat("\n-- meanme on model without me():\n")
r <- tryCatch(stancode(bf(y ~ x), data = d,
                       prior = set_prior("normal(0, 1)", class = "meanme")),
              error = function(e) conditionMessage(e))
print(substr(paste(r, collapse = " "), 1, 300))

cat("\n==================== L sample_prior = 'yes' with proper priors\n")
show("L yes on B", bf(y ~ me(x, sx) + me(z, sz)), prior = pr,
     sample_prior = "yes")
show("L2 yes on A with flat defaults", bf(y ~ me(x, sx)),
     sample_prior = "yes")
cat("\n-- 'only' on A with flat defaults:\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx)), data = d, sample_prior = "only"),
              error = function(e) conditionMessage(e))
print(substr(paste(r, collapse = " "), 1, 400))
cat("\n-- 'only' on A with proper priors:\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx)), data = d, sample_prior = "only",
                       prior = c(set_prior("normal(0, 2)", class = "meanme"),
                                 set_prior("exponential(1)",
                                           class = "sdme"))),
              error = function(e) conditionMessage(e))
cat(grep("me|target", strsplit(as.character(r), "\n")[[1]], value = TRUE),
    sep = "\n")
cat("\n-- 'only' on B with meanme/sdme proper and corme default:\n")
r <- tryCatch(stancode(bf(y ~ me(x, sx) + me(z, sz)), data = d,
                       sample_prior = "only",
                       prior = c(set_prior("normal(0, 2)", class = "meanme"),
                                 set_prior("exponential(1)",
                                           class = "sdme"))),
              error = function(e) conditionMessage(e))
cat(grep("me|target", strsplit(as.character(r), "\n")[[1]], value = TRUE),
    sep = "\n")
