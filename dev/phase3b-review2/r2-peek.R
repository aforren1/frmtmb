fs <- c(list.files("dev/phase3b-log/recov5", "^(cens|contfix)-", full.names=TRUE)[1:3], list.files("dev/phase3b-log/recov4","^cens-",full.names=TRUE)[1:3])
for (f in fs) {r<-readRDS(f); cat(f, r$elapsed, "\n")}
r<-readRDS(fs[4]); print(r$confint); print(names(r))
