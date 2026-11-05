source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
a <- readRDS(file.path(REV_OUT, "readers-lane.rds")); b <- readRDS(file.path(REV_OUT, "readers-base.rds"))
for (m in c("gp","grF","rep","aniso","cmcF","two","pois","by")) { cat("--",m,"\nlane "); print(unlist(a[[m]]$mfx)); cat("base "); print(unlist(b[[m]]$mfx)); cat("emtr lane/base "); print(c(unlist(a[[m]]$emtr[,2:3]), unlist(b[[m]]$emtr[,2:3]))); cat("mpred se ratio "); print(range(a[[m]]$mpred$std.error/b[[m]]$mpred$std.error)) }
