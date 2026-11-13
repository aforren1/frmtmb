.libPaths(c("C:/Users/adf44/source/r/wt-priortab-lib", "C:/Users/adf44/source/r/rellib-r20", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
F <- asNamespace("frmtmb")
for (lk in c("softplus", "squareplus")) {
  fam <- F$brmsfamily("gaussian", link_sigma = lk)
  sp <- F$frm(y ~ z, data = transform(data.frame(z = 1:20 / 10), y = z + sin(z)), family = fam, dry_run = "objective")
  li <- NULL
  for (lp in sp$frame$linpreds) if (identical(lp$dpar, "sigma")) li <- lp$link
  cat(lk, "\n")
  if (is.null(li)) { cat(" link object not found via spec$lp\n"); next }
  inv <- li$linkinv
  print(c(at_minus_inf = tryCatch(inv(-Inf), error = function(e) NaN), at_m700 = inv(-700), at_m1400 = inv(-1400)))
  print(F$link_image(li))
}
