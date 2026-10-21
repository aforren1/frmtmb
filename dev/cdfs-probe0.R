.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (f in c("pgamma","pbeta","pbinom","pnbinom","ppois","pt")) { cat("== RTMB::", f, "\n"); if (exists(f, asNamespace("RTMB"))) print(args(get(f, asNamespace("RTMB")))) }
for (f in c("pt","pskewnorm","pexgauss","pvm","pnbinom2","pzinbinom","pzipois","plaplace","pbeta2")) { cat("== RTMBdist::", f, "\n"); print(get(f, asNamespace("RTMBdist"))) }
