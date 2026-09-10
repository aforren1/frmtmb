# Probes that decide the design, run BEFORE any change.
#
# Seed 4242 throughout. Three questions:
#
#  1. What does rlddm() do today when it is handed ndt_group()?
#  2. Do the addition-term VALUES keep an attribute the coercion put on
#     them, all the way to the family and to newdata? That is the only
#     route by which the group LABEL could reach the floor lookup, so it
#     decides whether the label-keying defect filed against this item
#     can be closed at all.
#  3. Where does an aterm value get subset, and does the attribute
#     survive that?

source(file.path("C:/Users/adf44/source/r/frmtmb-wt-rlddm",
                 "dev/rlddm-scripts/rlddm-prelude.R"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})

set.seed(4242)

d <- frm_task_design("bandit2arm", n_subject = 4, n_trial = 40, seed = 1)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = 0.2, bias = 0.5), seed = 1)[[1L]]
cat("\nrows:", nrow(s), " min(rt):", format(min(s$rt), digits = 9), "\n")
cat("per-subject min(rt):\n")
print(tapply(s$rt, s$id, min))

# ---- 1. ndt_group() on rlddm today
f_g <- bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(id) ~ 1,
          drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)
cat("\n--- ndt_group() on rlddm(), before the change ---\n")
print(tryCatch(frm(f_g, family = rlddm(subject = id, trial = trial),
                   data = s, dry_run = "frame"),
               error = function(e) conditionMessage(e)))

# ---- 2. does an attribute put on by a coercion survive?
frmtmb_register_aterm("zzlab", arity = 1L, coerce = function(x) {
  structure(as.numeric(as.integer(factor(x))), zz_labels = as.character(x))
})
f_z <- bf(rt | dec(choice) + reward(pay1, pay2) + zzlab(id) ~ 1,
          drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)
fr <- frm(f_z, family = rlddm(subject = id, trial = trial), data = s,
          dry_run = "frame")
av <- fr$frame$aterm_values[["rt"]]
cat("\n--- attribute survival, in-sample frame ---\n")
cat("names(av):", paste(names(av), collapse = ", "), "\n")
cat("attributes(av$zzlab):",
    paste(names(attributes(av[["zzlab"]])), collapse = ", "), "\n")
cat("first 5 labels:",
    paste(utils::head(attr(av[["zzlab"]], "zz_labels"), 5),
          collapse = ", "), "\n")

# and on the newdata path, which is where predict() re-evaluates it
nd <- s[1:10, ]
anew <- tryCatch(
  frmtmb:::aterms_for_newdata(fr$spec$responses[["rt"]], nd),
  error = function(e) paste("ERR:", conditionMessage(e)))
cat("\n--- attribute survival, aterms_for_newdata ---\n")
if (is.list(anew)) {
  cat("names:", paste(names(anew), collapse = ", "), "\n")
  cat("attributes(zzlab):",
      paste(names(attributes(anew[["zzlab"]])), collapse = ", "), "\n")
} else {
  cat(anew, "\n")
}

# ---- 3. subsetting drops it, which is what the engine does per trial
v <- av[["zzlab"]]
cat("\n--- attribute after v[1:3] ---\n")
cat("attributes:", paste(names(attributes(v[1:3])), collapse = ", "),
    "| length ", length(attributes(v[1:3])), "\n")
