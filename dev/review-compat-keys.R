# Reviewer scratch (round 2): registered_aterm_of() / aterm_spelling() /
# aterm_base() over every registered term's keys, old rule against new,
# and constructed collisions in a scratch registry.
.libPaths(c("C:/Users/adf44/source/r/wt-compat-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (p in c("frmtmb.survival", "frmtmb.spatial", "frmtmb.latent",
            "frmtmb.sample", "frmtmb.learn", "frmtmb.eam",
            "frmtmb.coupling", "frmtmb.ode")) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
ns <- asNamespace("frmtmb")
reg_env <- ns$frmtmb_aterm_registry
old_of <- function(nm) {
  reg <- reg_env$reg
  if (!length(reg)) return(NULL)
  if (!is.null(reg[[nm]])) return(reg[[nm]])
  base <- sub("[0-9]+$", "", nm)
  e <- reg[[base]]
  if (!is.null(e) && e$arity > 1L) e else NULL
}
nm_of <- function(e) if (is.null(e)) "NULL" else e$name
keys_of <- function(e) {
  if (e$arity > 1L) paste0(e$name, seq_len(e$arity)) else e$name
}
show <- function(label) {
  cat("\n==", label, "\n")
  reg <- reg_env$reg
  cat(sprintf("%-12s arity %d option %s\n", names(reg),
              vapply(reg, `[[`, 1L, "arity"),
              vapply(reg, function(e) is.function(e$option), NA)), sep = "")
  keys <- unique(c(unlist(lapply(reg, keys_of)), "vint1", "vint12",
                   "vreal3", "trunc_lb", "cens_y2", "dec2", "reward3",
                   "payoff5", "stage2", "stage23", "stage"))
  for (k in keys) {
    cat(sprintf("%-12s new %-10s old %-10s base %-10s spell %s\n", k,
                nm_of(ns$registered_aterm_of(k)), nm_of(old_of(k)),
                ns$aterm_base(k), ns$aterm_spelling(k)))
  }
}
show("registry with all eight extensions")

# constructed collisions, in a scratch copy of the registry
saved <- reg_env$reg
put <- function(name, arity) {
  reg_env$reg[[name]] <- list(name = name, arity = as.integer(arity),
                              coerce = as.numeric, option = NULL)
}
reg_env$reg <- list()
put("a", 12); put("a1", 2)        # arity > 9 next to a digit-named term
put("pay", 2); put("payoff2", 3)  # a prefix of another
put("x", 3); put("x2", 1)         # arity-1 digit name inside x's keys
put("q1", 2); put("q", 1)         # digit name with arity > 1 next to base
show("constructed collisions (keys a11, a12, x2 are claimed twice)")
# what the parser writes for the colliding spellings
for (f in c("y | x(u, v, w) + x2(z) ~ 1", "y | x2(z) + x(u, v, w) ~ 1",
            "y | a1(u, v) + a(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12) ~ 1")) {
  r <- tryCatch(names(ns$parse_response(stats::as.formula(f))$aterms),
                error = function(e) paste("ERR:", conditionMessage(e)))
  cat(f, "->", paste(r, collapse = ","), "\n")
}
reg_env$reg <- saved
# does registration refuse a name whose keys another term already uses?
r <- tryCatch({frmtmb_register_aterm("reward1", arity = 1L); "registered"},
              error = function(e) paste("ERR:", conditionMessage(e)))
cat("\nregister reward1 beside reward (arity 2):", r, "\n")
r <- tryCatch({frmtmb_register_aterm("vint1", arity = 1L); "registered"},
              error = function(e) paste("ERR:", conditionMessage(e)))
cat("register vint1 beside core vint():", r, "\n")
reg_env$reg <- saved

t0 <- proc.time()[["elapsed"]]
for (i in 1:20) invisible(frm_compat("dec()", "thres()"))
cat("\nfrm_compat(a, b), mean of 20, all loaded:",
    round((proc.time()[["elapsed"]] - t0) / 20, 3), "s\n")
t0 <- proc.time()[["elapsed"]]
for (i in 1:5) invisible(ns$frmtmb_compat_rules_tbl())
cat("frmtmb_compat_rules_tbl(), mean of 5:",
    round((proc.time()[["elapsed"]] - t0) / 5, 3), "s\n")
