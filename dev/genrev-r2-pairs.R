# Every S3method(g, cls) on a borrowed name in frmtmb's own table must
# have S3method(owner::g, cls) for EACH owner, or the method goes
# unreachable the moment that owner loads.  Read from NAMESPACE; column
# 4 of parseNamespaceFile()'s matrix is the delayed-registration owner.
m <- parseNamespaceFile("frmtmb-wt-generics", "C:/Users/adf44/source/r")$S3methods
e <- new.env(); sys.source("R/generic-owners.R", envir = e, keep.source = FALSE)
own <- e$frm_generic_owners
local_rows <- which(is.na(m[, 4]) & m[, 1] %in% names(own))
miss <- character()
for (i in local_rows) for (o in own[[m[i, 1]]]) {
  if (!any(m[, 1] == m[i, 1] & m[, 2] == m[i, 2] & m[, 4] %in% o))
    miss <- c(miss, sprintf("%s.%s has no %s:: twin", m[i, 1], m[i, 2], o))
}
cat("borrowed-name registrations in frmtmb's own table:", length(local_rows), "\n")
cat("lacking an owner-table twin:", length(miss), "\n")
if (length(miss)) cat(paste0("  ", miss), sep = "\n")
