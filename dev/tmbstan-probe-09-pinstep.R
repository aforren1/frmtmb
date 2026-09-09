# lane tmbstan, probe 09: the rewritten pin step's own arithmetic and
# its three-state verdict, run against the runner's REAL repository
# table and against four model.hpp states.
#
# The repository pair is the one the runner prints in
# setup-r-dependencies under "Repo status", quoted in
# dev/reviews/2026-09-09-tmbstan.md:
#   1  RSPM  https://packagemanager.posit.co/cran/__linux__/noble/latest
#   2  CRAN  https://cran.rstudio.com
#
# The step's own source is extracted from the YAML rather than
# retyped, so this cannot drift from what will run.
yaml_path <- paste0("C:/Users/adf44/source/r/frmtmb-wt-tmbstan/",
                    ".github/workflows/check-frmtmb-sample.yaml")
src <- readLines(yaml_path)
i <- grep("^      - name: Pin a tmbstan", src)
j <- grep("^        shell: Rscript \\{0\\}", src)
j <- min(j[j > i])
body <- src[(i + 2L):(j - 1L)]
body <- sub("^          ", "", body)

# take only the two definitions and the repository arithmetic, which
# is everything before the install
cut <- grep("^install.packages", body)[1L]
head_src <- body[seq_len(cut - 1L)]

RUNNER <- c(RSPM = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
            CRAN = "https://cran.rstudio.com")

run_head <- function(repos, rspm_env = "") {
  e <- new.env(parent = globalenv())
  old_o <- options(repos = repos)
  old_e <- Sys.getenv("RSPM")
  Sys.setenv(RSPM = rspm_env)
  on.exit({
    options(old_o)
    if (nzchar(old_e)) Sys.setenv(RSPM = old_e) else Sys.unsetenv("RSPM")
  })
  out <- tryCatch({
    eval(parse(text = paste(head_src, collapse = "\n")), envir = e)
    list(ok = TRUE, pin = get("pin", envir = e), env = e)
  }, error = function(err) list(ok = FALSE, msg = conditionMessage(err)))
  out
}

cat("== B1: the repository the step reads ==\n")
r1 <- run_head(RUNNER)
cat("both entries present  -> ",
    if (r1$ok) paste("pin =", r1$pin) else paste("STOP:", r1$msg), "\n")

r2 <- run_head(RUNNER["CRAN"])
cat("CRAN only, no RSPM    -> ",
    if (r2$ok) paste("pin =", r2$pin) else
      paste("STOP:", substr(r2$msg, 1, 46), "..."), "\n")

r3 <- run_head(RUNNER["CRAN"], rspm_env = RUNNER[["RSPM"]])
cat("RSPM only in the env  -> ",
    if (r3$ok) paste("pin =", r3$pin) else paste("STOP:", r3$msg), "\n")

r4 <- run_head(c(RSPM = "https://packagemanager.posit.co/cran/__linux__/noble/2026-01-01",
                 CRAN = RUNNER[["CRAN"]]))
cat("RSPM already dated    -> ",
    if (r4$ok) paste("pin =", r4$pin) else
      paste("STOP:", substr(r4$msg, 1, 46), "..."), "\n")

r5 <- run_head(c(RSPM = paste0(RUNNER[["RSPM"]], "/")))
cat("RSPM trailing slash   -> ",
    if (r5$ok) paste("pin =", r5$pin) else paste("STOP:", r5$msg), "\n")

cat("\n== B2: the three-state verdict on four model.hpp states ==\n")
verdict <- get("verdict", envir = r1$env)
marker_line <-
  "    lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));"
clean_line <- "    lp_accum__.add(custom_func::custom_func(y));"
d <- tempfile("hppstates"); dir.create(d)
states <- list(
  absent = file.path(d, "does-not-exist.hpp"),
  `empty path` = "",
  empty = local({ p <- file.path(d, "e.hpp"); file.create(p); p }),
  patched = local({
    p <- file.path(d, "p.hpp"); writeLines(clean_line, p); p
  }),
  unpatched = local({
    p <- file.path(d, "u.hpp"); writeLines(c(clean_line, marker_line), p)
    p
  }))
for (nm in names(states)) {
  v <- verdict(states[[nm]])
  cat(sprintf("  %-11s -> verdict = %-9s bad = %s  job = %s\n", nm, v,
              !identical(v, "patched"),
              if (identical(v, "patched")) "continues" else "STOPS"))
}
unlink(d, recursive = TRUE)
