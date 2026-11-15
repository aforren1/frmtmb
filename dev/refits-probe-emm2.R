source("dev/refits-lib.R")
f <- deparse(emmeans::emmeans)
i <- grep("try[(]|ref_grid[(]|stop[(]", f)
for (k in i) cat(k, ":", f[k], "\n")
set.seed(6)
dj <- data.frame(x = rnorm(50), f = factor(rep(c("a", "b"), 25)))
dj$y <- 1 + dj$x + rnorm(50)
fj <- frm(y ~ x + f, data = dj)
show <- function(r) list(class = class(r), msg = conditionMessage(r))
go("ref_grid dpar unknown", show(tryCatch(
  emmeans::ref_grid(fj, dpar = "nope"), error = function(e) e)))
go("emmeans dpar unknown", show(tryCatch(
  emmeans::emmeans(fj, ~ f, dpar = "nope"), error = function(e) e)))
go("ref_grid data given", show(tryCatch(
  emmeans::ref_grid(fj, dpar = "nope", data = dj), error = function(e) e)))
go("ref_grid frmtmb_error handler", tryCatch(
  emmeans::ref_grid(fj, dpar = "nope"),
  frmtmb_error = function(e) paste("caught frmtmb_error:", e$message)))
go("try() text", {
  r <- try(emmeans::ref_grid(fj, dpar = "nope"), silent = TRUE)
  as.character(r)
})
