# Which objects of the base-priority packages that a session may not
# have loaded are NOT functions? base_r_function() decides those
# packages by exists() alone (forcing a promise would load the package),
# so these are the names it would call functions wrongly. Run in a
# throwaway process: forcing the promises here loads the packages.
for (p in c("splines", "stats4", "grid", "tools", "parallel", "compiler",
            "datasets", "tcltk")) {
  base <- file.path(find.package(p, lib.loc = .Library), "R", p)
  if (!file.exists(paste0(base, ".rdx"))) {
    cat(p, ": no lazy-load database\n"); next
  }
  e <- new.env(parent = emptyenv())
  lazyLoad(base, envir = e)
  nm <- ls(e, all.names = TRUE)
  nf <- nm[!vapply(nm, function(n) is.function(get(n, e)), NA)]
  cat(sprintf("%s: %d objects, %d not functions: %s\n", p, length(nm),
              length(nf), paste(nf, collapse = " ")))
}
