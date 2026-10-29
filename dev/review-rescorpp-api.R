## Review rescorpp: predict_rescor_draw() pairing of dps with rspecs, and
## predict() of a multivariate model with an ordinal response.
##   Rscript dev/review-rescorpp-api.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
source(file.path(wt, "dev/rescorpp-data.R"))
d <- rescorpp_data("g3")
fit <- quiet(frm(eval(rescorpp_frm_forms[["g3"]]), family = "gaussian", data = d))
dp <- eval_dpars(fit)
ok <- rep(TRUE, nrow(d))
set.seed(1)
a <- predict_rescor_draw(fit, fit$spec$responses[c("y3", "y1")], dp, ok)
cat("subset rspecs (y3, y1) with the full eval_dpars() list:\n")
cat("  names returned:", names(a), "\n")
cat("  mean of 'y3' draw", round(mean(a$y3), 3), " vs mu_y3 mean",
    round(mean(dp$y3$mu), 3), " vs mu_y1 mean", round(mean(dp$y1$mu), 3), "\n")
cat("  sd of 'y3' residual", round(sd(a$y3 - dp$y1$mu), 3), " sigma_y3",
    round(dp$y3$sigma[1], 3), " sigma_y1", round(dp$y1$sigma[1], 3), "\n")

dso <- readRDS(file.path(rv_out, "ds-ord-lane.rds"))
set.seed(2)
p_all <- quiet(predict(dso, ndraws = 200))
set.seed(2)
p_o1 <- quiet(predict(dso, resp = "o1", ndraws = 200))
cat("\nord: predict(ds) dim", dim(p_all), " layer names", dimnames(p_all)[[3]], "\n")
cat("predict(ds)[1:3, , 'o1']:\n"); print(p_all[1:3, , "o1"])
cat("predict(ds, resp = 'o1')[1:3, ]:\n"); print(p_o1[1:3, ])
show("predict(ds, resp = c('o1','y2')) dim", dim(predict(dso, resp = c("o1", "y2"), ndraws = 50)))
show("fitted(ds) ord", dim(fitted(dso, ndraws = 50)))
