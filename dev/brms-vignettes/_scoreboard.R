# Emit the Headline and Per-vignette sections of
# dev/brms-vignette-audit.md from the recorded CSVs.
#
# Why a generator: those sections are pure arithmetic over the bv() rows,
# and hand-typed counts went stale twice while the audit was being
# written. Run this and paste, or diff it against the file to check that
# the prose still matches the record.
#
#   BV_OUT=<dir> Rscript dev/brms-vignettes/_scoreboard.R

out <- Sys.getenv("BV_OUT", unset = tempdir())
d <- do.call(rbind, lapply(list.files(out, "[.]csv$", full.names = TRUE),
                           utils::read.csv, stringsAsFactors = FALSE))
d$edge[is.na(d$edge) | d$edge == ""] <- "CLEAN"
n_all <- nrow(d)
secs <- sum(d$secs)
d <- d[d$kind != "data", ]
d$path <- ifelse(grepl("^SAMPLE:", d$label), "SAMPLE", "ML")

# brms-coexistence is this audit's own instrument, not a shipped
# vignette, so it is scored apart from the headline.
inst <- "brms-coexistence"
v <- d[d$vignette != inst, ]
cls <- c("CLEAN", "SPELLING", "BEHAVIOR", "MISSING", "REFUSAL")

cat("## Headline\n\n")
cat("**", n_all, " translated calls across ",
    length(unique(d$vignette)), " scripts, ", round(secs),
    " seconds of compute.**\n\n", sep = "")

for (p in c("ML", "SAMPLE")) {
  s <- v[v$path == p, ]
  m <- s[s$kind == "model", ]
  po <- s[s$kind == "post", ]
  cat("### ", p, " path, ten vignettes\n\n", sep = "")
  cat("| | calls | port unchanged | share |\n|---|---|---|---|\n")
  f <- function(lbl, x) cat("| ", lbl, " | ", nrow(x), " | ",
                            sum(x$edge == "CLEAN"), " | ",
                            round(100 * sum(x$edge == "CLEAN") / nrow(x)),
                            "% |\n", sep = "")
  f("model calls", m)
  f("post-processing calls", po)
  f("**all calls**", rbind(m, po))
  cat("\n")
}

md <- v[v$kind == "model", ]
md$base <- sub("^(ML|SAMPLE): ", "", md$label)
u <- tapply(md$edge == "CLEAN", md$base, any)
cat("Distinct model calls across both paths: ", length(u),
    "; clean on at least one path: ", sum(u), ".\n\n", sep = "")

cat("## Per vignette\n\n")
cat("Model calls, then post-processing calls, as\n",
    "`total | clean spelling behavior missing refusal`.\n\n", sep = "")
for (p in c("ML", "SAMPLE")) {
  cat("### ", p, " path\n\n| vignette | model | post |\n|---|---|---|\n",
      sep = "")
  dp <- d[d$path == p, ]
  for (vv in c(sort(setdiff(unique(d$vignette), inst)), inst)) {
    s <- dp[dp$vignette == vv, ]
    cell <- function(x) {
      if (!nrow(x)) return("0")
      paste0(nrow(x), " \\| ",
             paste(vapply(cls, function(cc) sum(x$edge == cc), 1L),
                   collapse = " "))
    }
    nm <- if (vv == inst) paste0("*", vv, "*") else vv
    cat("| ", nm, " | ", cell(s[s$kind == "model", ]), " | ",
        cell(s[s$kind == "post", ]), " |\n", sep = "")
  }
  cat("\n")
}
