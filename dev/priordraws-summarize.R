# Counts and ranges for dev/priordraws-findings.md, generated from the
# logs in dev/priordraws-log/ rather than typed.
lg <- function(f) readLines(file.path("dev/priordraws-log", f))
num <- function(x) as.numeric(x)
cat("## names (names-frm.txt)\n")
cat(grep("^SUMMARY", lg("names-frm.txt"), value = TRUE), "\n")
p1 <- grep("KS D=", lg("evidence-1.txt"), value = TRUE)
pv <- num(sub(".* p=([0-9.]+).*", "\\1", p1))
models <- unique(sub("^ +([A-Z0-9]+) .*", "\\1", p1))
cat(sprintf("## part 1: %d columns of %d models, KS p from %.3f to %.3f; below 0.05: %d\n",
            length(pv), length(models), min(pv), max(pv), sum(pv < 0.05)))
p2 <- grep("KS D=", lg("evidence-2.txt"), value = TRUE)
pv2 <- num(sub(".* p=([0-9.]+)$", "\\1", p2))
cat(sprintf("## part 2: %d columns, two-sample KS p from %.3f to %.3f; below 0.05: %d\n",
            length(pv2), min(pv2), max(pv2), sum(pv2 < 0.05)))
p5 <- grep("KS D=", lg("evidence-5.txt"), value = TRUE)
pv5 <- num(sub(".* p=([0-9.]+)$", "\\1", p5))
cat(sprintf("## part 5: brms's own draws vs exact, %d columns, p from %.3f to %.3f\n",
            length(pv5), min(pv5), max(pv5)))
p3 <- grep("z mean", lg("evidence-3.txt"), value = TRUE)
nums <- function(s) {
  unlist(lapply(regmatches(s, gregexpr("-?[0-9]+[.][0-9]+", s)), num))
}
z3 <- nums(sub(".*z mean", "", p3))
cat(sprintf("## part 3: %d z-scores, max |z| %.2f\n", length(z3),
            max(abs(z3))))
p4 <- grep("brms ER", lg("evidence-4.txt"), value = TRUE)
z4 <- num(sub(".* z +(-?[0-9.]+)$", "\\1", p4))
cat(sprintf("## part 4: %d hypotheses, z from %.2f to %.2f\n", length(z4),
            min(z4), max(z4)))
p6 <- grep("^  seed", lg("evidence-6.txt"), value = TRUE)
for (side in c("brms", "frmtmb")) {
  s <- grep(paste0(" ", side, " "), p6, value = TRUE)
  z <- nums(sub(".*[|] z", "", sub("[|] ess_bulk.*", "", s)))
  cat(sprintf("## part 6: %s %d z-scores against the exact prior, max |z| %.2f\n",
              side, length(z), max(abs(z))))
}
idf <- list.files("dev/priordraws-log", "^identity-(base|lane)-[0-9]+[.]txt$")
if (length(idf)) {
  rows <- do.call(rbind, lapply(idf, function(f) {
    L <- lg(f)[-1]
    data.frame(file = f, model = sub(" .*", "", L),
               md5 = sub("^[^ ]+ ([^ ]+) .*", "\\1", L))
  }))
  for (m in unique(rows$model)) {
    r <- rows[rows$model == m, ]
    cat(sprintf("## identity %s: %d runs, distinct md5 %d (base %d, lane %d)\n",
                m, nrow(r), length(unique(r$md5)),
                length(unique(r$md5[grepl("base", r$file)])),
                length(unique(r$md5[grepl("lane", r$file)]))))
  }
}
