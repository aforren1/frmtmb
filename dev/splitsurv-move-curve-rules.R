# Move the curve-inference compat rules from frmtmb.spline's zzz.R into
# core's R/compat.R, after the frm_lp_basis() block. Run once.
rules <- readLines("dev/splitsurv-curve-rules.txt")
rules <- sub("^  ## ---- curve inference -+$",
             "  ## frm_curve() and its difference curve, which read frm_lp_basis()",
             rules)
rules <- sub("^  ## ---- the difference curve -+$", "", rules)
rules <- rules[nzchar(rules)]
fix <- c(
  "The case this package is for." = "The case frm_curve() is for.",
  "since this package reads frmtmb::frm_lp_basis()" =
    "since frm_curve() reads frm_lp_basis()",
  "the unit-perturbation rebuild this package used" =
    "the unit-perturbation rebuild frm_curve() used",
  "and this package used to call one" = "and frm_curve() used to call one",
  "the linearity probe this package used to run" =
    "the linearity probe frm_curve() used to run",
  "the joint precision this package inverts" =
    "the joint precision frm_curve() inverts")
for (k in names(fix)) {
  hit <- grepl(k, rules, fixed = TRUE)
  if (sum(hit) != 1L) stop("expected one hit for: ", k)
  rules <- sub(k, fix[[k]], rules, fixed = TRUE)
}
stopifnot(!any(grepl("this package", rules)))
f <- "R/compat.R"
x <- readLines(f)
i <- grep("^  r[(]\"frm_lp_basis\", \"predict\", \"works\",", x)
stopifnot(length(i) == 1L)
x <- c(x[1:(i + 1L)], "", rules, x[(i + 2L):length(x)])
j <- grep("^             \"frm_lp_basis\"[)], f,$", x)
stopifnot(length(j) == 1L)
x[j] <- "             \"frm_lp_basis\", \"frm_curve\", \"frm_curve_contrast\"), f,"
writeLines(x, f)
cat("inserted", length(rules), "lines\n")
