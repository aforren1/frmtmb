# genrev round 2: can the active binding's memo be POISONED?
#
# frm_bind_generic() memoizes (owner namespace, owner generic) and
# revalidates only on namespace IDENTITY.  If the binding is read while
# an owner's namespace is registered but its exports are not yet set,
# getExportedValue() fails, NULL is memoized against that namespace,
# and nothing invalidates it after the load completes.
#
# loadNamespace() runs packageEvent(pkg, "onLoad") hooks beside .onLoad
# (line 367) and namespaceExport() only at line 520, so any hook that
# reads frmtmb's binding for a name that owner exports reaches exactly
# that window.  A user hook is the construction; the question is what
# a read in that window leaves behind.
LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
say <- function(...) cat(sprintf(...))
suppressMessages(library(frmtmb))
own <- function(f) environmentName(topenv(environment(f)))
seen <- NULL
setHook(packageEvent("loo", "onLoad"), function(...) {
  seen <<- c(loaded = isNamespaceLoaded("loo"),
             exported = tryCatch({ getExportedValue("loo", "loo"); TRUE },
                                 error = function(e) FALSE),
             binding_owner = own(get("loo", envir = asNamespace("frmtmb"))))
})
suppressMessages(loadNamespace("brms"))
say("LIB %s\n", LIB)
say("inside the loo onLoad hook: loaded=%s exported=%s frmtmb loo -> %s\n",
    seen[["loaded"]], seen[["exported"]], seen[["binding_owner"]])
f <- get("loo", envir = globalenv())
tb <- get(".__S3MethodsTable__.", envir = environment(f))
say("after load: frmtmb loo -> %s ; loo.brmsfit reachable %s\n",
    own(f), exists("loo.brmsfit", envir = tb, inherits = FALSE))
say("waic (never read in the window) -> %s\n",
    own(get("waic", envir = globalenv())))
cat("GENREVDONE\n")
