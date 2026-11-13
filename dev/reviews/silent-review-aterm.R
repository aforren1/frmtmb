# Review of item 4675: extension load orders, reload, and edge keys.
libs <- c("C:/Users/adf44/source/r/wt-silent-lib",
          "C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
ld <- function(p) tryCatch({
  suppressPackageStartupMessages(loadNamespace(p)); "ok"
}, error = function(e) paste("ERROR", conditionMessage(e)))
ord <- c("frmtmb.survival", "frmtmb.learn", "frmtmb.eam", "frmtmb.latent",
         "frmtmb.coupling", "frmtmb.spatial")
for (p in ord) cat("load", p, ld(p), "\n")
for (p in rev(ord)) {
  r <- tryCatch({unloadNamespace(p); "unloaded"},
                error = function(e) paste("unload ERROR", conditionMessage(e)))
  cat(p, r, "\n")
}
for (p in rev(ord)) cat("reload", p, ld(p), "\n")
for (p in ord) cat("again", p, ld(p), "\n")
reg <- frmtmb:::frmtmb_aterm_registry$reg
cat("registered:", paste(vapply(reg, function(e) paste0(e$name, "/",
                                                        e$arity), ""),
                         collapse = " "), "\n")
# register the same name twice with a different arity (a replacement)
r1 <- tryCatch({frmtmb_register_aterm("zz", arity = 2L); "ok"},
               error = function(e) conditionMessage(e))
r2 <- tryCatch({frmtmb_register_aterm("zz", arity = 3L); "ok"},
               error = function(e) conditionMessage(e))
cat("zz/2", r1, " zz/3 again", r2, "\n")
for (nm in c("zz1", "zz3", "zz4", "vint", "vreal", "mi_s", "trunc_l",
             "thres_g", "se_s", "cens_y")) {
  r <- tryCatch({frmtmb_register_aterm(nm); "registered"},
                error = function(e) substr(conditionMessage(e), 1, 120))
  cat("register", nm, ":", r, "\n")
}
