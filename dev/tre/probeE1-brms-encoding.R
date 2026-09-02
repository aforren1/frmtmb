# Probe E. What exactly does brms's gr(dist = "student") mean, and what
# does it refuse? The spelling is not being invented here, so the
# encoding is the specification, and every combination frmtmb might have
# to guard has to be asked of brms first.
#
# Run: Rscript dev/tre/probeE1-brms-encoding.R
sink("dev/tre/probeE1.txt")
suppressMessages(library(brms))
set.seed(1)
d <- data.frame(y = rnorm(120), x = rnorm(120), z = rnorm(120),
                g = rep(1:20, each = 6), h = rep(1:6, 20),
                by = rep(c("a", "b"), each = 60))
A <- diag(20); dimnames(A) <- list(1:20, 1:20)

cat("brms ", as.character(packageVersion("brms")), "\n\n", sep = "")

cat("gr() formals:\n")
print(names(formals(brms::gr)))
cat("\ndist is matched against: ")
print(eval(formals(brms::gr)$dist))
cat("\n")

show <- function(label, form, ...) {
  cat("\n---- ", label, " ----\n", sep = "")
  r <- tryCatch(make_stancode(form, data = d, ...),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (startsWith(r[1], "ERROR:")) { cat(r, "\n"); return(invisible()) }
  ln <- strsplit(r, "\n")[[1]]
  keep <- grep(paste0("df_|udf_|dfm_|inv_chi_square|r_1|r_2|sd_1|sd_2|",
                      "L_1|z_1|Lcov|student"), ln, value = TRUE)
  cat(paste(unique(keep), collapse = "\n"), "\n")
}

show("scalar intercept", y ~ x + (1 | gr(g, dist = "student")))
show("correlated slopes", y ~ x + (1 + z | gr(g, dist = "student")))
show("uncorrelated (cor = FALSE)",
     y ~ x + (1 + z | gr(g, dist = "student", cor = FALSE)))
show("with cov = A", y ~ x + (1 | gr(g, dist = "student", cov = A)),
     data2 = list(A = A))
show("with by =", y ~ x + (1 | gr(g, dist = "student", by = by)))
show("two terms sharing id =",
     y ~ x + (1 | q | gr(g, dist = "student")) +
       (0 + z | q | gr(g, dist = "student")))
show("student on one group, gaussian on another",
     y ~ x + (1 | gr(g, dist = "student")) + (1 | h))
show("nested inside mm()",
     y ~ x + (1 | mm(g, h, dist = "student")))

cat("\n---- default prior on df ----\n")
print(default_prior(y ~ x + (1 | gr(g, dist = "student")), data = d))

cat("\n---- can df be fixed by a constant prior? ----\n")
r <- tryCatch({
  p <- prior(constant(3), class = "df", group = "g")
  sc <- make_stancode(y ~ x + (1 | gr(g, dist = "student")), data = d,
                      prior = p)
  ln <- strsplit(sc, "\n")[[1]]
  paste(grep("df_1", ln, value = TRUE), collapse = "\n")
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, "\n")

cat("\n---- standata: is dist visible in the DATA block? ----\n")
sd1 <- make_standata(y ~ x + (1 | gr(g, dist = "student")), data = d)
sd0 <- make_standata(y ~ x + (1 | g), data = d)
cat("student names: ", paste(names(sd1), collapse = " "), "\n")
cat("gaussian names: ", paste(names(sd0), collapse = " "), "\n")
cat("identical: ", identical(lapply(sd1, as.vector),
                             lapply(sd0, as.vector)), "\n")
cat("=> the choice is encoded in the STANCODE only; make_standata\n")
cat("   cannot be used as the structural cross-check here.\n")

cat("\n---- how brms names the parameter in the fit summary ----\n")
b <- brms::brmsterms(y ~ x + (1 + z | gr(g, dist = "student")))
print(b$dpars$mu$re)

sink()
cat("done\n")
