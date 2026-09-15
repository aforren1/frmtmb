source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")

# Attribution, not a guess: dev/generics-adoptparts.R times the pieces.
# isNamespaceLoaded 0.20 us, asNamespace 0.92, identical 0.27,
# getExportedValue 0.37, and the whole body inline 1.60 us. The
# measured binding cost was 8.89 us, so about 7 us was not in the body
# at all: it is the `tryCatch()` wrapped around `asNamespace()`,
# which installs a handler stack on every access. `asNamespace()`
# cannot fail behind an `isNamespaceLoaded()` that just returned TRUE,
# so the guard comes off the hot path and stays on the cold one, where
# `getExportedValue()` runs once per namespace instance.
sub1("R/generic-owners.R",
paste0("    if (is.null(own)) return(fallback)\n",
       "    ons <- tryCatch(asNamespace(own), error = function(e) NULL)\n",
       "    if (is.null(ons)) return(fallback)\n"),
paste0("    if (is.null(own)) return(fallback)\n",
       "    # no tryCatch on this line: it runs on every access, a\n",
       "    # handler stack costs about 7 us here, and asNamespace()\n",
       "    # cannot fail behind the isNamespaceLoaded() above\n",
       "    ons <- asNamespace(own)\n"))
cat("DONE\n")
