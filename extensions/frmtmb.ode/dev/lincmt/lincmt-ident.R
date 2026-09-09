# frm_lincmt() against frm_ode() over the schedule space.
#
# The solver is run at atol = rtol = 1e-12 rather than at its default
# 1e-8, so that what is measured is the closed form and not the
# integrator. Both a scale-relative and a pointwise-relative error are
# reported, because the two disagree exactly where the solver's own
# absolute tolerance dominates.
source("lincmt-src.R")

mk_dyn <- function(ncmt, depot) {
  if (!depot && ncmt == 1L) return(function(t, y, p) list(c(-p[1] * y[1])))
  if (!depot && ncmt == 2L) {
    return(function(t, y, p) list(c(
      -(p[1] + p[2]) * y[1] + p[3] * y[2],
      p[2] * y[1] - p[3] * y[2])))
  }
  if (!depot && ncmt == 3L) {
    return(function(t, y, p) list(c(
      -(p[1] + p[2] + p[4]) * y[1] + p[3] * y[2] + p[5] * y[3],
      p[2] * y[1] - p[3] * y[2],
      p[4] * y[1] - p[5] * y[3])))
  }
  if (ncmt == 1L) {
    return(function(t, y, p) list(c(
      -p[6] * y[1], p[6] * y[1] - p[1] * y[2])))
  }
  if (ncmt == 2L) {
    return(function(t, y, p) list(c(
      -p[6] * y[1],
      p[6] * y[1] - (p[1] + p[2]) * y[2] + p[3] * y[3],
      p[2] * y[2] - p[3] * y[3])))
  }
  function(t, y, p) list(c(
    -p[6] * y[1],
    p[6] * y[1] - (p[1] + p[2] + p[4]) * y[2] + p[3] * y[3] + p[5] * y[4],
    p[2] * y[2] - p[3] * y[3],
    p[4] * y[2] - p[5] * y[4]))
}

st_names <- function(ncmt, depot) {
  c(if (depot) "depot", "central",
    if (ncmt >= 2L) "peripheral1", if (ncmt == 3L) "peripheral2")
}

pk_list <- function(ncmt, depot, p) {
  z <- list(ke = p[["ke"]])
  if (ncmt >= 2L) { z$k12 <- p[["k12"]]; z$k21 <- p[["k21"]] }
  if (ncmt == 3L) { z$k13 <- p[["k13"]]; z$k31 <- p[["k31"]] }
  if (depot) z$ka <- p[["ka"]]
  z
}

pvec <- function(p) c(p[["ke"]], p[["k12"]], p[["k21"]], p[["k13"]],
                      p[["k31"]], p[["ka"]])

# One case: returns the two error measures and the trajectory scale.
run_case <- function(ncmt, depot, p, times, ev = NULL, init = NULL,
                     t0 = 0, group = NULL, escale = 1, n_ss = 20L,
                     out = "central", tol = 1e-12) {
  nst <- st_names(ncmt, depot)
  ic <- if (depot) 2L else 1L
  i0 <- rep(list(0), ncmt + as.integer(depot))
  if (!is.null(init)) {
    if (!is.null(init$depot)) i0[[1L]] <- init$depot
    if (!is.null(init$central)) i0[[ic]] <- init$central
  }
  a <- frm_lincmt(parms = pk_list(ncmt, depot, p), times = times,
                  group = group, ncmt = ncmt, depot = depot,
                  output = out, t0 = t0, init = init, events = ev,
                  event_scale = escale, n_ss = n_ss)
  b <- frm_ode(mk_dyn(ncmt, depot), init = i0, times = times,
               parms = as.list(pvec(p)), group = group, t0 = t0,
               states = nst, output = out, events = ev,
               event_scale = escale, n_ss = as.integer(n_ss),
               atol = tol, rtol = tol)
  s <- max(abs(b))
  list(scale = s,
       rel_scale = if (s > 0) max(abs(a - b)) / s else max(abs(a - b)),
       rel_point = max(abs(a - b) / pmax(abs(b), 1e-300)),
       min_over_scale = if (s > 0) min(abs(b)) / s else NA_real_)
}

P0 <- list(ke = 0.2, k12 = 0.4, k21 = 0.1, k13 = 0.05, k31 = 0.01,
           ka = 1.1)
tt <- c(0, 0.5, 1, 2, 4, 8, 12, 18, 24, 30, 36, 48)

cases <- list()
add <- function(name, ...) cases[[length(cases) + 1L]] <<-
  c(list(name = name), list(...))

for (nc in 1:3) for (dp in c(FALSE, TRUE)) {
  tag <- paste0(nc, "cmt", if (dp) "+depot" else "")
  cmt_in <- if (dp) "depot" else "central"
  add(paste(tag, "init bolus"), ncmt = nc, depot = dp, p = P0,
      times = tt, init = if (dp) list(depot = 100) else
        list(central = 100))
  add(paste(tag, "single dose"), ncmt = nc, depot = dp, p = P0,
      times = tt,
      ev = data.frame(time = 2, state = cmt_in, value = 100))
  add(paste(tag, "addl/ii"), ncmt = nc, depot = dp, p = P0, times = tt,
      ev = data.frame(time = 0, state = cmt_in, value = 100, ii = 8,
                      addl = 5L))
  add(paste(tag, "ss"), ncmt = nc, depot = dp, p = P0, times = tt,
      ev = data.frame(time = 0, state = cmt_in, value = 100, ii = 8,
                      ss = TRUE))
  add(paste(tag, "ss + addl"), ncmt = nc, depot = dp, p = P0, times = tt,
      ev = data.frame(time = c(0, 8), state = cmt_in, value = 100,
                      ii = c(8, 8), addl = c(0L, 5L),
                      ss = c(TRUE, FALSE)))
  add(paste(tag, "infusion single"), ncmt = nc, depot = dp, p = P0,
      times = tt,
      ev = data.frame(time = 1, state = "central", value = 100,
                      duration = 3))
  add(paste(tag, "infusion addl"), ncmt = nc, depot = dp, p = P0,
      times = tt,
      ev = data.frame(time = 0, state = "central", value = 100,
                      duration = 3, ii = 8, addl = 5L))
  add(paste(tag, "infusion ss"), ncmt = nc, depot = dp, p = P0,
      times = tt,
      ev = data.frame(time = 0, state = "central", value = 100,
                      duration = 3, ii = 8, ss = TRUE))
  add(paste(tag, "reset"), ncmt = nc, depot = dp, p = P0, times = tt,
      init = if (dp) list(depot = 100) else list(central = 100),
      ev = data.frame(time = c(12, 12), state = c(NA, cmt_in),
                      value = c(0, 100),
                      method = c("reset", "add")))
  add(paste(tag, "t0 = 3"), ncmt = nc, depot = dp, p = P0,
      times = tt[tt >= 3], t0 = 3,
      ev = data.frame(time = 4, state = cmt_in, value = 100))
  if (dp) {
    add(paste(tag, "depot output"), ncmt = nc, depot = dp, p = P0,
        times = tt, out = "depot",
        ev = data.frame(time = 0, state = "depot", value = 100, ii = 8,
                        addl = 5L))
  }
}

# degenerate rate constants
deg <- list(
  "ka == ke" = modifyList(P0, list(ka = 0.2)),
  "ka == ke exactly, both 1" = modifyList(P0, list(ka = 1, ke = 1)),
  "ka - ke = 1e-9" = modifyList(P0, list(ka = 0.2 + 1e-9)),
  "ka - ke = 1e-13" = modifyList(P0, list(ka = 0.2 + 1e-13)),
  "k12 -> 0, k21 == ke" = modifyList(P0, list(k12 = 1e-10, k21 = 0.2)),
  "k12 -> 0, k21 == ke == ka" =
    modifyList(P0, list(k12 = 1e-10, k21 = 0.2, ka = 0.2)),
  "k31 == k21" = modifyList(P0, list(k31 = 0.1, k13 = 0.4)),
  "all peripheral rates equal" =
    modifyList(P0, list(k12 = 0.3, k21 = 0.3, k13 = 0.3, k31 = 0.3)),
  "everything equal" =
    modifyList(P0, list(ke = 0.3, k12 = 0.3, k21 = 0.3, k13 = 0.3,
                        k31 = 0.3, ka = 0.3)))
for (nm in names(deg)) {
  for (nc in 1:3) {
    add(paste0(nc, "cmt+depot ", nm), ncmt = nc, depot = TRUE,
        p = deg[[nm]], times = tt,
        ev = data.frame(time = 0, state = "depot", value = 100, ii = 8,
                        addl = 3L))
    add(paste0(nc, "cmt+depot ", nm, ", ss"), ncmt = nc, depot = TRUE,
        p = deg[[nm]], times = tt,
        ev = data.frame(time = 0, state = "depot", value = 100, ii = 8,
                        ss = TRUE))
  }
}

# groups, a per-group schedule and an estimated dose scale
gd <- data.frame(id = factor(rep(1:3, each = 6)),
                 time = rep(c(0, 1, 4, 8, 12, 24), 3))
add("3 groups, per-group schedule, event_scale", ncmt = 2, depot = TRUE,
    p = P0, times = gd$time, group = gd$id, escale = c(0.6, 0.9, 1.2)[
      as.integer(gd$id)],
    ev = data.frame(group = as.character(rep(1:3, each = 2)),
                    time = rep(c(0, 8), 3), state = "depot",
                    value = rep(c(100, 50, 200), each = 2),
                    ii = 8, addl = rep(c(1L, 2L, 0L), each = 2)))

res <- data.frame()
for (cs in cases) {
  nm <- cs$name
  cs$name <- NULL
  r <- tryCatch(do.call(run_case, cs), error = function(e) e)
  if (inherits(r, "error")) {
    res <- rbind(res, data.frame(case = nm, scale = NA, rel_scale = NA,
                                 rel_point = NA, min_rel = NA,
                                 note = conditionMessage(r)))
  } else {
    res <- rbind(res, data.frame(case = nm, scale = r$scale,
                                 rel_scale = r$rel_scale,
                                 rel_point = r$rel_point,
                                 min_rel = r$min_over_scale, note = ""))
  }
}
res$note <- NULL
# grouped by SCHEDULE KIND, which is the case name with its model
# shape stripped: that is what the findings document reports
res$kind <- sub("^[123]cmt([+]depot)? ", "", res$case)
g <- aggregate(cbind(rel_scale, rel_point) ~ kind, data = res, FUN = max)
n <- aggregate(rel_scale ~ kind, data = res, FUN = length)
g$n <- n$rel_scale[match(g$kind, n$kind)]
cat("--- worst per schedule kind ---
")
print(g[order(-g$rel_scale), c("kind", "n", "rel_scale", "rel_point")],
      right = FALSE, row.names = FALSE, digits = 3)
cat("--- top 12 by scale-relative error ---
")
print(head(res[order(-res$rel_scale), ], 12), right = FALSE,
      row.names = FALSE, digits = 3)
cat("--- top 12 by pointwise relative error ---
")
print(head(res[order(-res$rel_point), ], 12), right = FALSE,
      row.names = FALSE, digits = 3)
cat("\nWORST rel_scale:", format(max(res$rel_scale, na.rm = TRUE)),
    " WORST rel_point:", format(max(res$rel_point, na.rm = TRUE)), "\n")
cat("cases:", nrow(res), " errors:", sum(nzchar(res$note)), "\n")
