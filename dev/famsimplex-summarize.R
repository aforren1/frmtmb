# Generate the numbered blocks of dev/famsimplex-findings.md from the
# lane's logs and result files, so no count in the record is typed.
# Output: dev/famsimplex-summary.md (pasted into the findings verbatim).
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
out <- character(0)
emit <- function(...) out <<- c(out, paste0(...))
block <- function(title, lines) {
  emit("<!-- BEGIN generated: ", title, " -->")
  emit("```")
  for (l in lines) emit(l)
  emit("```")
  emit("<!-- END generated -->")
  emit("")
}
rd <- function(f) readLines(file.path(wt, f), warn = FALSE)

# 1. Stan log densities
s <- rd("famsimplex-stan-lpdf.log")
i <- grep("WORST per family", s)
block("Stan lpdf grid (dev/famsimplex-stan-lpdf.R, seed 20261013)",
      c(s[i:(i + 4L)], grep("^rows compared", s, value = TRUE)))

# 2. brms_lp_check rows
l <- rd("famsimplex-lpcheck.log")
lp <- grep("^LPCHECK", l, value = TRUE)
rows <- grep("^[0-9] +row", l, value = TRUE)
v <- as.numeric(sub("^LPCHECK const (\\S+) .*$", "\\1", lp))
g <- as.numeric(sub("^.* grad (\\S+) .*$", "\\1", lp))
block("brms_lp_check (dev/famsimplex-lpcheck.R)",
      c(paste(sub("^[0-9] +", "", trimws(sub(" +[0-9]+ +[0-9]+ +FALSE$",
                                               "", rows))),
              "|", lp),
        sprintf("rows %d; max |const| %.3g; max |grad| %.3g", length(lp),
                max(abs(v)), max(g))))

# 3. agreement with brms fits
a <- rd("famsimplex-agree.log")
i <- grep("^SUMMARY", a)
block("ML vs brms posterior (dev/famsimplex-agree.R)",
      c(a[i:length(a)]))

# 4. recovery
for (cfg in c("dir", "dirre", "ln", "lnre", "dir2")) {
  r <- rd(paste0("famsimplex-recovery-", cfg, ".log"))
  block(paste0("recovery ", cfg, " (dev/famsimplex-recovery.R ", cfg,
               " 1001 100)"), r[-1L])
}

# 5. brms tolerance experiment
t2 <- rd("famsimplex-brms-tol2.log")
keep <- grep("^== row|^draws:|no more scalars|sampling not done", t2,
             value = TRUE)
block("brms dirichlet, one row off the simplex (dev/famsimplex-brms-tol2.R)",
      trimws(keep))
t1 <- rd("famsimplex-brms-tol.log")
block("brms, row 1 off by 1e-7 (dev/famsimplex-brms-tol.R)",
      trimws(grep("^row 1|^== |lncor|sampled|no more scalars", t1,
                  value = TRUE)))

# 6. edge cases against brms standata
e <- rd("famsimplex-edge.log")
keep <- grep("^=== |^brms |^ERROR|^frm fits|^MESSAGE", e, value = TRUE)
block("response checks beside brms standata() (dev/famsimplex-edge.R)",
      substr(keep, 1, 150))

# 7. pins on rellib-r8
for (arm in c("before", "after")) {
  block(paste0("pins, ", arm, " (dev/famsimplex-pins.R)"),
        substr(rd(paste0("famsimplex-pins-", arm, ".log")), 1, 150))
}
block("mixture mu1 ~ x (dev/famsimplex-mixmu1.R)",
      c(rd("famsimplex-mixmu1-before.log"), rd("famsimplex-mixmu1-after.log")))

# 8. test runs
for (tag in c("before1", "before2", "new3", "new4")) {
  f <- file.path(wt, "famsimplex-tests", paste0(tag, ".log"))
  if (file.exists(f)) {
    block(paste0("test run ", tag, " (dev/famsimplex-run.sh)"),
          sub(" [|] lib: .*; frmtmb [0-9.]+ ", " | core ",
              readLines(f, warn = FALSE)))
  }
}
writeLines(out, file.path(wt, "famsimplex-summary.md"))
cat("wrote", length(out), "lines\n")

# 9. the whole-suite run, per package, and the reruns
out <- character(0)
f <- file.path(wt, "famsimplex-tests", "all1.log")
r <- readLines(f, warn = FALSE)
jobs <- readLines(file.path(wt, "famsimplex-tests-all.jobs"))
pk <- sub(" .*$", "", jobs)
res <- grep("^RESULT", r, value = TRUE)
num <- function(k) {
  as.integer(sub(paste0("^.* ", k, "=([0-9]+).*$"), "\\1", res))
}
tab <- aggregate(cbind(pass = num("pass"), fail = num("fail"),
                       err = num("err"), skip = num("skip"),
                       warn = num("warn")) ~ pk, FUN = sum)
tab$files <- as.integer(table(pk)[tab$pk])
libs <- sum(grepl("wt-famsimplex-lib", res))
block("whole suite, gated (dev/famsimplex-run.sh all1)",
      c(utils::capture.output(print(tab, row.names = FALSE)),
        tail(r, 1), sprintf("RESULT lines naming the lane library: %d", libs),
        grep("fail=[1-9]|err=[1-9]", res, value = TRUE) |>
          sub(pattern = " [|] lib: .*$", replacement = "")))
for (tag in c("extra1", "ob326")) {
  block(paste0("test run ", tag, " (dev/famsimplex-run.sh / -ob.sh)"),
        sub(" [|] lib: .*; frmtmb [0-9.]+ ", " | core ",
            readLines(file.path(wt, "famsimplex-tests",
                                paste0(tag, ".log")), warn = FALSE)))
}
pf <- rd("famsimplex-postfit-agree.log")
block("fitted() and predict() against brms (dev/famsimplex-postfit-agree.R)",
      pf[-(1:2)])
writeLines(out, file.path(wt, "famsimplex-summary2.md"))
cat("wrote summary2", length(out), "lines\n")
