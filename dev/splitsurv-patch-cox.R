f <- "extensions/frmtmb.survival/R/cox-family.R"
x <- readLines(f)
sub1 <- function(old, new) {
  i <- which(x == old)
  if (length(i) != 1L) stop("not once: ", old)
  x[i] <<- new
}
sub1("  # item 4). `chart` is the response's `bhaz_chart` addition term, which",
     "  # item 4). `chart` is the response's `sbhaz_chart` addition term, which")
sub1("  # frame (mo_chart_frame()).",
     "  # frame (frmtmb::mo_chart_frame(), for simplex_pars below).")
sub1("    links = list(mu = mu_link(link, \"cox\")),",
     "    links = list(mu = link),")
x <- gsub("aterms[[\"bhaz_chart\"]]", "aterms[[\"sbhaz_chart\"]]", x,
          fixed = TRUE)
i <- which(x == "      \"distribution\"")
stopifnot(length(i) == 1L, x[i + 1L] == "    )", x[i + 2L] == "  )")
j <- which(x == "  fam[[\"aterm_data\"]] <- function(y, aterms) {")
k <- which(x == "  fam[[\"cox_sbhaz\"]] <- sbhaz")
body <- x[(j + 1L):(k - 2L)]
stopifnot(x[k - 1L] == "  }")
new_tail <- c(
  "    ),",
  "    # the baseline bases are data, not parameters: they are built once",
  "    # from the observed times and ride with the response's addition",
  "    # terms",
  "    aterm_data = function(y, aterms) {",
  paste0("  ", body),
  "    },",
  "    # the baseline weights are a simplex, which core searches, charts",
  "    # and priors as brms's class sbhaz",
  "    simplex_pars = c(sbhaz = \"sbhaz_raw\")",
  "  )")
x <- c(x[1:i], new_tail, x[k:length(x)])
writeLines(x, f)
