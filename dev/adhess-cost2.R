# Lane adhess: the cost of the changed component, arms interleaved in
# one process (lane library). Per model (the models of
# dev/adhess-cost.R, same seeds), at the fit's own point:
#   old   fd_hessian(at, gr) at its default steps: 0.71.0's outer Hessian
#   new   the lane's: obj$he(at) symmetrized without random effects
#         (uncached, the worst case), fd_hessian() with
#         outer_fd_steps() with them
#   ctrl  `old` again, which must report 1.0 against `old`
# Each block repeats its call until it has run 1.2 s; the time per call
# is the minimum over 5 rounds, the rounds interleaved A B C. Also
# printed: whether the fit-time code had already cached the AD Hessian
# at that point (then the lane's covariance costs no new Hessian).
#   Rscript dev/adhess-cost2.R <model>
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
model <- commandArgs(TRUE)[1]
src <- parse("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-cost.R")
# the model list `mk` of dev/adhess-cost.R, evaluated as written there
mk_expr <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                    identical(e[[2]], as.name("mk")), as.list(src))
eval(mk_expr[[1]])
cs <- mk[[model]]()
fit <- suppressWarnings(frm(cs$f, family = cs$fam, data = cs$d))
obj <- fit$obj
at <- ns$sdr_outer_point(obj)
hit <- obj$env$.frmtmb_he
cached <- !is.null(hit) && identical(hit$p, unname(as.numeric(at)))
re <- length(obj$env$random) > 0
st <- ns$obj_state_save(obj)
old <- function() ns$fd_hessian(at, obj$gr)$H
new <- if (!re) function() {
  H <- obj$he(at)
  (H + t(H)) / 2
} else {
  eps <- ns$outer_fd_steps(fit$frame, names(fit$opt$par))
  function() ns$fd_hessian(at, obj$gr, eps = eps)$H
}
block <- function(f) {
  k <- 0L
  t0 <- proc.time()[["elapsed"]]
  repeat {
    f()
    k <- k + 1L
    el <- proc.time()[["elapsed"]] - t0
    if (el >= 1.2) return(el / k)
  }
}
tt <- matrix(NA, 5, 3, dimnames = list(NULL, c("old", "new", "ctrl")))
for (r in 1:5) {
  tt[r, "old"] <- block(old)
  tt[r, "new"] <- block(new)
  tt[r, "ctrl"] <- block(old)
}
ns$obj_state_restore(obj, st)
m <- apply(tt, 2, min)
cat(sprintf(paste0("%-18s RE %-5s np %3d n %6d | cached at fit %-5s | ",
                   "old %.4f s new %.4f s ctrl %.4f s | new/old %.3f ",
                   "ctrl/old %.3f\n"),
            model, re, length(at), nrow(cs$d), cached, m[["old"]],
            m[["new"]], m[["ctrl"]], m[["new"]] / m[["old"]],
            m[["ctrl"]] / m[["old"]]))
