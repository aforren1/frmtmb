# Short check of fix round 1, M5 and m1: the floor arms on the repaired warning.
#   Rscript dev/review-gpgap-floorf1.R
L <- unlist(lapply(c("dev/review-gpgap-out/floorf1-rep.txt",
                     "dev/review-gpgap-out/floorf1-seeds.txt"), readLines, warn = FALSE))
L <- L[startsWith(L, "FLOOR")]
g <- function(p) sub(p, "\\1", L)
tb <- data.frame(
  tag = g("^FLOOR (\\S+) .*"),
  ratio = as.numeric(g(".* ratio (\\S+) \\| brms.*")),
  rel = as.numeric(g(".*sigma/brms (\\S+), .*")),
  dse = suppressWarnings(as.numeric(g(".*\\|dlog\\|/SE (\\S+), .*"))),
  sdb = as.numeric(g(".*brms: sigma \\S+ sdgp (\\S+) rho.*")),
  sdf = as.numeric(g(".*frmnug: sigma \\S+ sdgp (\\S+) rho.*")),
  gap = as.numeric(g(".*brms: sigma \\S+ sdgp \\S+ rho \\S+ ll (\\S+) \\| frmnug.*")) -
    as.numeric(g(".*frmnug: sigma \\S+ sdgp \\S+ rho \\S+ ll (\\S+) \\| sigma/brms.*")),
  w = grepl("nugwarn TRUE", L),
  stringsAsFactors = FALSE)
print(tb[tb$tag != "seeds", ], row.names = FALSE)
s <- tb[tb$tag == "seeds", ]
cat(sprintf("seeds rows %d, warned %d, warned below 1 SE %d, not warned at or above 1 SE %d\n",
            nrow(s), sum(s$w), sum(s$w & s$dse < 1, na.rm = TRUE),
            sum(!s$w & s$dse >= 1, na.rm = TRUE)))
print(aggregate(cbind(w, dse) ~ ratio, s, mean), digits = 3)
