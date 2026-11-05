# Fill the generated numbers of the short check and append it to the review.
s <- readLines("C:/Users/adf44/AppData/Local/Temp/2/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/predread-shortcheck.md", warn = FALSE)
s <- paste(s, collapse = "\n")
O <- "dev/review-predread-out/r1/"
tot <- function(d) {
  fs <- list.files(paste0(O, d), full.names = TRUE)
  r <- unlist(lapply(fs, function(f) grep("^RESULT", readLines(f, warn = FALSE), value = TRUE)))
  g <- function(k) sum(as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", r)))
  c(files = length(r), of = length(fs), pass = g("pass"), fail = g("fail"), err = g("err"), skip = g("skip"), warn = g("warn"))
}
a <- tot("merge-tests")
b <- tot("merge-ext")
res <- sprintf(paste0("**Result.** Core and frmtmb.sample, one file per process, gated: %d of %d files, %d pass, %d fail, %d errors, %d skips, %d escaped warnings. ",
                      "`test-latres-brms.R` 278 pass, `test-nl-brms.R` 91 pass, both to the end, fail 0. The seven other suites (on the 23:55 snapshot, which differs from this one in two nlfix files of core): %d of %d files, %d pass, %d fail, %d errors. ",
                      "This is a green suite, not a clean merge: R1, S4c and the items below pass every test. `R CMD check` was not run on the merge, and nlfix has moved since this snapshot only if its files are newer than 00:32."),
               a["files"], a["of"], a["pass"], a["fail"], a["err"], a["skip"], a["warn"], b["files"], b["of"], b["pass"], b["fail"], b["err"])
s <- sub("MERGE_RESULT", res, s, fixed = TRUE)
s <- sub("Snapshot at 00:3x with", "Snapshot at 00:32 with", s, fixed = TRUE)
s <- sub("SESIM_RESULT", "the lane's numbers reproduce. Coverage of the mean at the drawn latent at `sdx = 1`: 95.5 and 98.0 percent. Negative lower ends: 0. Truth / estimate 0.993 to 1.002.", s, fixed = TRUE)
s <- sub("by up to 3.5 times in width", "by up to 3.3 times in width", s, fixed = TRUE)
s <- sub("Bernoulli 0.423, 0.574, 0.706, 0.642 against 0.421, 0.572, 0.713, 0.647.",
         "Bernoulli `fitted` 0.423, 0.574, 0.714, 0.647 against draws 0.421, 0.572, 0.713, 0.647.", s, fixed = TRUE)
stopifnot(!grepl("MERGE_RESULT|SESIM_RESULT", s))
f <- "dev/reviews/2026-11-02-predread.md"
cat(s, "\n", sep = "", file = f, append = TRUE)
print(a)
print(b)
