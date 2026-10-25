# Fix round 1: the core tests for F1 (the map check, a single sigma),
# F2 (class priors on a bounded parameter) and F14 (a seam function of
# `...`) in tests/testthat/test-autocor-seam.R. Rscript <this> <root>
root <- commandArgs(TRUE)[1]
p <- file.path(root, "tests/testthat/test-autocor-seam.R")
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * as.vector(aux$L %*% u)
    })
}",
"    # u is a vector of innovations or a matrix of them, one column per
    # set; %*% keeps the shape
    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * (aux$L %*% u)
    })
}")
rep1("test_that(\"sigma_varies = FALSE refuses a predicted sigma\", {
  local_frmtmb_registries()
  register_myfix(sigma_varies = FALSE)
  sd_ <- seam_data()",
"test_that(\"sigma_varies = FALSE refuses a predicted sigma and passes one
           sigma\", {
  local_frmtmb_registries()
  # the map is handed one sigma, as ?frmtmb_register_autocor says, and
  # refuses anything else
  tt <- myfix_term()
  one <- function(f) {
    function(u, mu, sigma, theta, aux) {
      if (length(sigma) != 1L) stop(\"sigma has \", length(sigma), \" values\")
      f(u, mu, sigma, theta, aux)
    }
  }
  frmtmb_register_autocor(\"myfix\", args = \"V\", build = tt$build,
                          whiten = one(tt$whiten),
                          unwhiten = one(tt$unwhiten), sigma_varies = FALSE)
  sd_ <- seam_data()
  f1 <- frm(y ~ x + myfix(V), data = sd_$d, data2 = list(V = sd_$V))
  f0 <- frm(y ~ x + fcor(V), data = sd_$d, data2 = list(V = sd_$V))
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  expect_identical(fitted(f1), fitted(f0))
  expect_identical(simulate(f1, nsim = 2, seed = 4),
                   simulate(f0, nsim = 2, seed = 4))
  expect_identical(residuals(f1, type = \"pearson\"),
                   residuals(f0, type = \"pearson\"))")
rep1("  expect_s3_class(frm(y ~ x + myfix(V), data = sd_$d,
                      data2 = list(V = sd_$V)), \"frmtmb_fit\")
})",
"  expect_s3_class(frm(y ~ x + myfix(V), data = sd_$d,
                      data2 = list(V = sd_$V)), \"frmtmb_fit\")
})

test_that(\"a map whose two halves do not agree is refused at the frame\", {
  local_frmtmb_registries()
  tt <- myfix_term()
  reg <- function(name, whiten = tt$whiten, unwhiten = tt$unwhiten) {
    frmtmb_register_autocor(name, args = \"V\", build = tt$build,
                            whiten = whiten, unwhiten = unwhiten)
  }
  # whiten() that gives one innovation fewer than the rows
  reg(\"shortu\", whiten = function(y, mu, sigma, theta, aux) {
    w <- tt$whiten(y, mu, sigma, theta, aux)
    w$u <- w$u[-1L]
    w
  })
  # unwhiten() that is not the inverse of whiten(): the transpose of L
  reg(\"wrongun\", unwhiten = function(u, mu, sigma, theta, aux) {
    mu + sigma * (t(aux$L) %*% u)
  })
  # unwhiten() that flattens a matrix of innovations
  reg(\"flatun\", unwhiten = function(u, mu, sigma, theta, aux) {
    mu + sigma * as.vector(aux$L %*% u)
  })
  # whiten() whose log-Jacobian is one value per row
  reg(\"vecjac\", whiten = function(y, mu, sigma, theta, aux) {
    w <- tt$whiten(y, mu, sigma, theta, aux)
    w$logjac <- rep(w$logjac / length(y), length(y))
    w
  })
  sd_ <- seam_data()
  fit_with <- function(nm) {
    frm(stats::as.formula(paste0(\"y ~ x + \", nm, \"(V)\")), data = sd_$d,
        data2 = list(V = sd_$V))
  }
  expect_error(fit_with(\"shortu\"),
               \"'shortu' term's map is wrong: whiten() gave 23 innovation\",
               fixed = TRUE, class = \"frmtmb_error\")
  expect_error(fit_with(\"wrongun\"),
               \"'wrongun' term's map is wrong: whiten(unwhiten(u)) is not u\",
               fixed = TRUE, class = \"frmtmb_error\")
  expect_error(fit_with(\"flatun\"),
               \"unwhiten() of a matrix of innovations must give\",
               fixed = TRUE, class = \"frmtmb_error\")
  expect_error(fit_with(\"vecjac\"), \"must give one finite log-Jacobian\",
               fixed = TRUE, class = \"frmtmb_error\")
  # the absent case: core's fcor() and a correct map pass the same check
  # (frmtmb.spatial's sar() passes it in that package's test-sar.R)
  reg(\"rightmap\")
  expect_s3_class(fit_with(\"rightmap\"), \"frmtmb_fit\")
  expect_s3_class(fit_with(\"fcor\"), \"frmtmb_fit\")
})

test_that(\"a bounded parameter takes a prior by its name\", {
  local_frmtmb_registries()
  tt <- myfix_term()
  half <- function(n) seq_len(n) <= n / 2
  # the scale of the first half of the rows, on (0.5, 2)
  sc <- function(theta) 0.5 + 1.5 / (1 + exp(-theta[1]))
  frmtmb_register_autocor(
    \"bscale\", args = \"V\", build = tt$build, npar = 1L,
    natural = function(theta, aux) c(scale1 = sc(theta)),
    bounds = function(aux) cbind(lower = 0.5, upper = 2),
    whiten = function(y, mu, sigma, theta, aux) {
      tt$whiten(y, mu, sigma * sc(theta)^half(length(y)), theta, aux)
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      tt$unwhiten(u, mu, sigma * sc(theta)^half(NROW(u)), theta, aux)
    })
  sd_ <- seam_data()
  d <- sd_$d
  V <- sd_$V
  fit <- frm(y ~ x + bscale(V), data = d, data2 = list(V = V))
  nat1 <- function(f) {
    unname(frmtmb:::autocor_natural(f$estimates$thetaac,
                                    f$frame$autocor[[1L]]))
  }
  ap <- frmtmb_autocor_pars(fit)
  expect_identical(ap$name, \"scale1\")
  expect_identical(c(ap$lower, ap$upper), c(0.5, 2))
  dp <- default_prior(y ~ x + bscale(V), data = d, data2 = list(V = V))
  expect_true(\"scale1\" %in% dp$class)
  # a tight prior on scale1 pulls it to 1.8; the same prior written on
  # the internal parameter, with the log-Jacobian of the map, is the
  # same model
  pr <- frm(y ~ x + bscale(V), data = d, data2 = list(V = V),
            prior = set_prior(\"normal(1.8, 0.01)\", class = \"scale1\"))
  s1 <- nat1(pr)
  expect_lt(abs(s1 - 1.8) / 1.8, 0.01)
  # lb and ub narrow the interval
  pb <- frm(y ~ x + bscale(V), data = d, data2 = list(V = V),
            prior = set_prior(\"\", class = \"scale1\", lb = 1.5, ub = 1.9))
  expect_gte(nat1(pb), 1.5)
  # a parameter whose term states no interval is refused by name, with
  # the internal spelling to use instead
  frmtmb_register_autocor(
    \"uscale\", args = \"V\", build = tt$build, npar = 1L,
    natural = function(theta, aux) c(scale2 = exp(theta[1])),
    whiten = function(y, mu, sigma, theta, aux) {
      tt$whiten(y, mu, sigma * exp(theta[1] * half(length(y))), theta, aux)
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      tt$unwhiten(u, mu, sigma * exp(theta[1] * half(NROW(u))), theta, aux)
    })
  expect_error(frm(y ~ x + uscale(V), data = d, data2 = list(V = V),
                   prior = set_prior(\"normal(1, 1)\", class = \"scale2\")),
               \"whose term declares no interval\", class = \"frmtmb_error\")
  # bounds() that do not match natural() are refused at the frame
  frmtmb_register_autocor(
    \"badbound\", args = \"V\", build = tt$build, npar = 1L,
    natural = function(theta, aux) c(scale3 = exp(theta[1])),
    bounds = function(aux) cbind(lower = 0.5, upper = 2),
    whiten = function(y, mu, sigma, theta, aux) {
      tt$whiten(y, mu, sigma * exp(theta[1] * half(length(y))), theta, aux)
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      tt$unwhiten(u, mu, sigma * exp(theta[1] * half(NROW(u))), theta, aux)
    })
  expect_error(frm(y ~ x + badbound(V), data = d, data2 = list(V = V)),
               \"natural() is not the bounded map of bounds()\", fixed = TRUE,
               class = \"frmtmb_error\")
})")
rep1("    unwhiten = function(u, mu, sigma, theta, aux) {
      n <- length(u)
      tt$unwhiten(u, mu, sigma * exp(theta[1] * half(n)), theta, aux)
    })",
"    unwhiten = function(u, mu, sigma, theta, aux) {
      n <- NROW(u)
      tt$unwhiten(u, mu, sigma * exp(theta[1] * half(n)), theta, aux)
    })")
rep1("  ex(reg(\"ok1\", start = NA_real_, npar = 0L), \"start\")",
"  ex(reg(\"ok1\", start = NA_real_, npar = 0L), \"start\")
  ex(reg(\"ok1\", bounds = 2), \"bounds\")")
rep1("  reg(\"ok1\")
  expect_true(\"ok1\" %in% frmtmb:::ac_entry_names())",
"  reg(\"ok1\")
  expect_true(\"ok1\" %in% frmtmb:::ac_entry_names())
  # a slot that is a function of `...` takes any arguments
  frmtmb_register_autocor(\"ok2\", args = \"V\", build = tt$build,
                          whiten = function(...) tt$whiten(...),
                          unwhiten = function(...) tt$unwhiten(...))
  expect_true(\"ok2\" %in% frmtmb:::ac_entry_names())")
rep1("  expect_false(any(c(\"myfix\", \"twoscale\", \"ok1\") %in%
                     frmtmb:::ac_entry_names()))",
"  expect_false(any(c(\"myfix\", \"twoscale\", \"ok1\", \"ok2\", \"bscale\",
                     \"shortu\") %in% frmtmb:::ac_entry_names()))")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("edited", p, "\n")
