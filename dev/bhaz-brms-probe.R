# What brms 2.23.0 does with its bhaz() addition term and with
# cox(bhaz = ). No compile: standata(), stancode(), default_prior().
# Then compare brms's Zbhaz and Zcbhaz with frmtmb.survival's basis on
# the same times. Run with rellib-r12 (frmtmb.survival 0.1.0) behind
# the user library:
#   Rscript dev/bhaz-brms-probe.R > dev/bhaz-log/brms-probe.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), " splines2",
    format(packageVersion("splines2")), "\n")
set.seed(11); n <- 200; x <- rnorm(n)
g <- factor(sample(c("a", "b", "c"), n, TRUE))
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = g)
cat("range(time)", format(range(dd$time), digits = 17), "\n")

show <- function(lab, fo, family = cox()) {
  cat("\n=====", lab, "=====\n")
  r <- tryCatch({
    sd <- standata(fo, data = dd, family = family)
    sc <- stancode(fo, data = dd, family = family)
    dp <- default_prior(fo, data = dd, family = family)
    list(sd = sd, sc = sc, dp = dp)
  }, error = function(e) e, warning = function(w) w)
  if (inherits(r, "condition")) {
    cat(class(r)[1L], ":", conditionMessage(r), "\n")
    return(invisible(NULL))
  }
  sd <- r$sd
  cat("standata names with bhaz:",
      grep("bhaz", names(sd), value = TRUE), "\n")
  cat("Kbhaz", sd$Kbhaz, " dim Zbhaz", dim(sd$Zbhaz), "\n")
  if (!is.null(sd$ngrbhaz)) {
    cat("ngrbhaz", sd$ngrbhaz, " table(Jgrbhaz)",
        table(sd$Jgrbhaz), "\n")
    cat("con_sbhaz dim", dim(sd$con_sbhaz), "\n")
  }
  at <- attributes(sd$Zbhaz)
  for (a in c("degree", "knots", "Boundary.knots", "intercept",
              "derivs", "integral", "periodic")) {
    if (!is.null(at[[a]])) cat(" attr", a, ":",
                               format(at[[a]], digits = 17), "\n")
  }
  cat("stancode lines with bhaz:\n")
  cat(grep("bhaz", strsplit(r$sc, "\n")[[1]], value = TRUE), sep = "\n")
  cat("default_prior:\n")
  print(as.data.frame(r$dp)[, c("prior", "class", "coef", "group",
                                "source")])
  invisible(sd)
}

s0 <- show("cox(), no bhaz()", bf(time | cens(cens) ~ x))
s6 <- show("bhaz(df = 6)", bf(time | cens(cens) + bhaz(df = 6) ~ x))
kn <- unname(quantile(dd$time, c(0.3, 0.6)))
bk <- c(0, max(dd$time) * 1.1)
cat("\nknots used:", format(kn, digits = 17), " Boundary.knots used:",
    format(bk, digits = 17), "\n")
skb <- show("bhaz(knots =, Boundary.knots =)",
            bf(time | cens(cens) + bhaz(knots = kn, Boundary.knots = bk)
               ~ x))
sk <- show("bhaz(knots =) alone",
           bf(time | cens(cens) + bhaz(knots = kn) ~ x))
skd <- show("bhaz(df = 8, knots =): df and knots disagree",
            bf(time | cens(cens) + bhaz(df = 8, knots = kn) ~ x))
sb <- show("bhaz(Boundary.knots =) alone",
           bf(time | cens(cens) + bhaz(Boundary.knots = bk) ~ x))
sdi <- show("bhaz(df = 5, degree = 2, intercept = FALSE)",
            bf(time | cens(cens) + bhaz(df = 5, degree = 2,
                                        intercept = FALSE) ~ x))
bin <- c(0.5, max(dd$time) * 0.8)
sbin <- show("bhaz(Boundary.knots =) inside the data range",
             bf(time | cens(cens) + bhaz(Boundary.knots = bin) ~ x))
sg <- show("bhaz(gr = g)", bf(time | cens(cens) + bhaz(gr = g) ~ x))
sgd <- show("bhaz(gr = g, df = 4)",
            bf(time | cens(cens) + bhaz(gr = g, df = 4) ~ x))
show("cox(bhaz = list(df = 6))", bf(time | cens(cens) ~ x),
     family = tryCatch(brms::cox(bhaz = list(df = 6)),
                       error = function(e) {
                         cat("\nbrms::cox(bhaz = list(df = 6)) errors:",
                             conditionMessage(e), "\n")
                         cox()
                       }))
cat("args of brms::cox:", names(formals(brms::cox)), "\n")
show("gaussian with bhaz(df = 6)",
     bf(time | bhaz(df = 6) ~ x), family = gaussian())
show("bhaz() with no arguments", bf(time | cens(cens) + bhaz() ~ x))
show("bhaz(6) positional", bf(time | cens(cens) + bhaz(6) ~ x))

# ---- the same basis in frmtmb.survival 0.1.0 (rellib-r12) ----
ns <- asNamespace("frmtmb.survival")
cmp <- function(lab, sd, df, degree = 3L, intercept = TRUE) {
  sp <- ns$bhaz_spec(dd$time, df, degree, intercept)
  M <- ns$mspline_design(dd$time, sp)
  I <- ns$ispline_design(dd$time, sp)
  bz <- unclass(sd$Zbhaz); attributes(bz) <- list(dim = dim(sd$Zbhaz))
  cz <- unclass(sd$Zcbhaz); attributes(cz) <- list(dim = dim(sd$Zcbhaz))
  at <- attributes(sd$Zbhaz)
  cat(sprintf(paste0("%-40s dims %s/%s  max|dM|/max|M| %.3e  ",
                     "max|dI|/max|I| %.3e  knots equal %s  boundary ",
                     "equal %s\n"),
              lab, paste(dim(M), collapse = "x"),
              paste(dim(bz), collapse = "x"),
              max(abs(M - bz)) / max(abs(bz)),
              max(abs(I - cz)) / max(abs(cz)),
              isTRUE(all.equal(unname(sp$internal),
                               unname(as.numeric(at$knots)),
                               tolerance = 0)),
              isTRUE(all.equal(sp$boundary,
                               as.numeric(at$Boundary.knots),
                               tolerance = 0))))
}
cat("\n===== frmtmb.survival 0.1.0 basis against brms's =====\n")
cmp("cox() vs brms default", s0, 5L)
cmp("cox(df = 6) vs bhaz(df = 6)", s6, 6L)
cmp("cox(df = 5, degree = 2, intercept = FALSE)", sdi, 5L, 2L, FALSE)
