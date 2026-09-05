# esicar lane: dump brms's generated Stan code for the four car types.
lib_es <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
.libPaths(c(lib_es, "C:/Users/adf44/AppData/Local/R/win-library/4.6", .libPaths()))
suppressMessages(library(brms))
set.seed(19)
g <- expand.grid(r = 1:4, c = 1:4)
loc <- paste0("l", seq_len(nrow(g)))
W <- matrix(0, nrow(g), nrow(g), dimnames = list(loc, loc))
for (i in seq_len(nrow(g))) for (j in seq_len(nrow(g))) {
  if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
}
d <- data.frame(loc = factor(rep(loc, 3), levels = loc),
                x = rnorm(48), y = rnorm(48))
out <- "C:/Users/adf44/source/r/frmtmb-wt-esicar/dev/es-stancode"
dir.create(out, showWarnings = FALSE)
for (ty in c("escar", "esicar", "icar", "bym2")) {
  f <- bf(as.formula(paste0("y ~ x + car(W, gr = loc, type = \"", ty, "\")")))
  code <- brms::make_stancode(f, data = d, family = gaussian(),
                              data2 = list(W = W))
  writeLines(as.character(code), file.path(out, paste0(ty, ".stan")))
}
sdat <- brms::make_standata(bf(y ~ x + car(W, gr = loc, type = "esicar")),
                            data = d, family = gaussian(),
                            data2 = list(W = W))
cat("standata names:", paste(names(sdat), collapse = ", "), "\n")
cat("Nloc:", sdat$Nloc, " Nedges:", sdat$Nedges, "\n")
cat("brms version:", as.character(utils::packageVersion("brms")), "\n")
