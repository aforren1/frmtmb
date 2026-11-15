# Short check of lane gpfield's fix round, F1: gp_end_warn() said once on
# every route; the cache flag across saveRDS(), update() and refit().
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
ns <- asNamespace("frmtmb")
cat("has gp_end_warn:", exists("gp_end_warn", ns), "\n")
cnt <- function(w) {
  cl <- sub(":.*", "", w)
  tb <- table(cl[grepl("^frmtmb_gp", cl)])
  if (!length(tb)) return("none")
  paste(names(tb), tb, sep = " x", collapse = ", ")
}
fits <- list()
for (tg in c("sweep:mm_V:bernoulli:40:1", "sweep:exp_a:bernoulli:120:2",
             "gate:approx_a", "ng:student:step:80:1")) {
  r <- gf_row(tg)
  q <- gf_fit(r$fml, data = r$data, family = r$family)
  f <- q$fit
  fits[[tg]] <- list(fit = f, r = r)
  cat(sprintf("F1 %-28s conv %s | gp warnings: %s | all: [%s] | flag %s\n",
              tg, format(f$opt$convergence), cnt(q$warn),
              paste(sub(":.*", "", q$warn), collapse = ";"),
              format(isTRUE(f$cache$gp_end_said))))
}
# saveRDS / readRDS, then update() and refit() from the read fit
tg <- "sweep:mm_V:bernoulli:40:1"
f <- fits[[tg]]$fit
r <- fits[[tg]]$r
tf <- tempfile(fileext = ".rds")
saveRDS(f, tf)
f2 <- readRDS(tf)
cat("read fit flag:", isTRUE(f2$cache$gp_end_said), "size MB",
    round(file.size(tf) / 2^20, 2), "\n")
w <- character(0)
f3 <- withCallingHandlers(update(f2, newdata = r$data),
                          warning = function(e) {
                            w <<- c(w, paste0(class(e)[1], ": ",
                                              conditionMessage(e)))
                            invokeRestart("muffleWarning")
                          }, message = function(m) invokeRestart("muffleMessage"))
cat("F1 update(readRDS fit) | gp warnings:", cnt(w), "| ll",
    format(as.numeric(logLik(f3)), digits = 10), "\n")
w <- character(0)
f4 <- withCallingHandlers(refit(f2, newresp = r$data$y),
                          warning = function(e) {
                            w <<- c(w, paste0(class(e)[1], ": ",
                                              conditionMessage(e)))
                            invokeRestart("muffleWarning")
                          }, message = function(m) invokeRestart("muffleMessage"))
cat("F1 refit(readRDS fit, same y) | gp warnings:", cnt(w), "| ll",
    format(as.numeric(logLik(f4)), digits = 10), "| messages on it now:",
    length(ns$gp_end_messages(f4)), "\n")
# a non-nl fit: the fit-end block alone
d <- fits[["ng:student:step:80:1"]]$r$data
w <- character(0)
f5 <- withCallingHandlers(update(readRDS({saveRDS(fits[["ng:student:step:80:1"]]$fit, tf); tf})),
                          warning = function(e) {
                            w <<- c(w, paste0(class(e)[1], ": ",
                                              conditionMessage(e)))
                            invokeRestart("muffleWarning")
                          }, message = function(m) invokeRestart("muffleMessage"))
cat("F1 update(readRDS st80) | gp warnings:", cnt(w), "\n")
# gp_end_warn() twice on one fit: once
w <- character(0)
withCallingHandlers({ns$gp_end_warn(f2); ns$gp_end_warn(f2)},
                    warning = function(e) {
                      w <<- c(w, class(e)[1]); invokeRestart("muffleWarning")
                    })
cat("F1 gp_end_warn twice on read fit:", length(w), "warnings\n")
