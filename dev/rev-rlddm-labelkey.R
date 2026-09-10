# REVIEW, item 1.0b, attack 5: is the label really unreachable from
# frmtmb.eam, or only unreachable through the coercion's RETURN VALUE?
#
#   Rscript dev/rev-rlddm-labelkey.R <lib>
#
# Seed 4242. The lane's own construction is reproduced first: a
# coercion that attaches the labels as an attribute, read back in the
# fitted frame and on the newdata path. Then the case the lane did not
# construct: the coercion cannot RETURN the label, but it is CALLED
# with it, on both paths. A side map written there and read at lookup
# time needs no change to frmtmb at all, and that is what decides
# whether "not closable inside frmtmb.eam" is a fact or a preference.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("lib:", lib, " eam", format(packageVersion("frmtmb.eam")), "\n\n")

set.seed(4242)
d <- ddm_simulate(200, mu = 1.2, bs = 1.5, ndt = 0.25)
d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))

# A side map that lives where a package's own namespace environment
# would. Nothing about it needs frmtmb to change.
seen <- new.env(parent = emptyenv())
calls <- new.env(parent = emptyenv())
calls$n <- 0L
calls$where <- character(0)

frmtmb_register_aterm("ndt_group", arity = 1L, coerce = function(x) {
  lab <- as.character(x)
  co <- ndt_bound_key(x)
  calls$n <- calls$n + 1L
  calls$where <- c(calls$where, paste(sort(unique(lab)), collapse = "+"))
  for (i in seq_along(lab)) assign(format(co[i], digits = 17), lab[i],
                                   envir = seen)
  structure(co, zz_labels = lab)
})

fit <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
              bias = 0.5), family = wiener(), data = d)
av <- fit$frame$aterm_values[["rt"]]
att <- function(v) if (length(attributes(v))) {
  paste(names(attributes(v)), collapse = ", ")
} else "NONE"

cat("-- the lane's measurement, reproduced --\n")
cat("aterm names in the fitted frame :",
    paste(names(av), collapse = ", "), "\n")
cat("attributes on the fitted value  :", att(av[["ndt_group"]]), "\n")
nd <- d[1:6, ]
anew <- frmtmb:::aterms_for_newdata(frmtmb::single_response(fit), nd)
cat("attributes on the newdata value :", att(anew[["ndt_group"]]), "\n")
cat("attributes after v[1:3]         :", att(av[["ndt_group"]][1:3]), "\n")

cat("\n-- what the coercion itself SAW, on both paths --\n")
cat("coercion calls                  :", calls$n, "\n")
cat("labels each call was handed     :",
    paste(calls$where, collapse = " | "), "\n")
cat("side map size after both paths  :", length(ls(seen)), "\n")
cat("side map contents               :",
    paste(vapply(ls(seen), function(k) paste0(k, "->", get(k, seen)), ""),
          collapse = ", "), "\n")

# the residual, and whether the side map closes it: a label the fit
# never saw, whose CODE collides with one it did
cat("\n-- the residual, with and without the side map --\n")
code_a <- ndt_bound_key("a")
cat("code for 'a'                    :", format(code_a, digits = 17), "\n")
cat("a colliding label is detected by the side map:",
    !identical(get(format(code_a, digits = 17), seen), "zzz"), "\n")
cat("(that is, the map holds 'a' at that code, so a newdata row whose\n",
    " label is not 'a' but whose code is can be refused by comparing\n",
    " the label the coercion was just handed against the stored one)\n")

# and the cost of statefulness, stated rather than hidden: a SECOND fit
# on different labels shares the same map
d2 <- d
d2$g <- factor(rep(c("x", "y"), length.out = nrow(d2)))
fit2 <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
               bias = 0.5), family = wiener(), data = d2)
cat("\nafter a second, unrelated fit   :", length(ls(seen)),
    "entries in one process-wide map\n")
