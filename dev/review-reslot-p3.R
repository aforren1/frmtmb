# Reviewer probe 3: one slot per parameter, on structures the lane did
# not list. For each model: every position cov_class_coefs() gives a
# class must (1) be reached by that class's own spelling with the
# fields the refusal names, (2) be refused under class theta naming
# that class; every class theta row default_prior() lists must resolve,
# and every theta / thetar position must be in exactly one of the two.
# Usage: Rscript dev/review-reslot-p3.R after. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
set.seed(1)
G <- 8L; n <- 12L
dd <- data.frame(g = factor(rep(seq_len(G), each = n)),
                 t = rep(seq_len(n), G))
dd$h <- factor(rep(c("a", "b"), length.out = G))[dd$g]
dd$g2 <- factor(sample(letters[1:4], nrow(dd), TRUE))
dd$g3 <- factor(sample(letters[5:8], nrow(dd), TRUE))
dd$x <- rnorm(nrow(dd)); dd$z <- rnorm(nrow(dd))
dd$tf <- factor((dd$t - 1) %% 3)
dd$tn <- as.numeric(dd$t)
dd$y <- 1 + 0.5 * dd$x + rnorm(G, 0, 0.7)[dd$g] + rnorm(nrow(dd))
dd$y2 <- 0.5 * dd$y + rnorm(nrow(dd))
dd$y3 <- 0.3 * dd$y + rnorm(nrow(dd))
dd$yn <- 2 * exp(0.3 * dd$x) + rnorm(nrow(dd), 0, 0.3)
dd$ord <- cut(dd$y, 4, labels = FALSE)
dd$cat <- factor(cut(dd$y, 3, labels = c("A", "B", "C")))
A <- diag(G) * 0.5 + 0.5; dimnames(A) <- list(levels(dd$g), levels(dd$g))
d2 <- list(A = A)

msg_of <- function(expr) tryCatch({ force(expr); "" },
                                  error = function(e) conditionMessage(e))

reached <- function(des, pl) {
  r <- frmtmb:::resolve_priorlist(des, pl)
  pt <- des$frame[["par_template"]]
  nm <- unlist(lapply(r$entries, function(e) {
    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
  }))
  sort(unique(c(nm, names(r$lower), names(r$upper))))
}

dens <- c(sd = "normal(0, 1)", sdgp = "normal(0, 1)",
          lscale = "inv_gamma(2, 1)", cor = "lkj(2)", rescor = "lkj(2)")

probe <- function(label, f, family = NULL, data2 = list()) {
  cat("== ", label, "\n", sep = "")
  des <- tryCatch(frmtmb:::prior_design(f, dd, family, data2),
                  error = function(e) e)
  if (inherits(des, "error")) {
    cat("  DESIGN ERROR:", conditionMessage(des), "\n"); return(invisible())
  }
  cc <- frmtmb:::cov_class_coefs(des$spec, des$frame)
  tab <- tryCatch(default_prior(f, data = dd, family = family,
                                data2 = data2), error = function(e) e)
  if (inherits(tab, "error")) {
    cat("  TABLE ERROR:", conditionMessage(tab), "\n"); return(invisible())
  }
  pt <- des$frame[["par_template"]]
  all_th <- c(frmtmb:::par_template_names(pt[["theta"]] %||% numeric(0),
                                          "theta"),
              frmtmb:::par_template_names(pt[["thetar"]] %||% numeric(0),
                                          "thetar"))
  th_rows <- tab$coef[tab$class == "theta" & nzchar(tab$coef)]
  cat("  theta positions:", length(all_th), " with a class:", nrow(cc),
      " theta rows:", length(th_rows), " classwide theta row:",
      any(tab$class == "theta" & !nzchar(tab$coef)), "\n")
  both <- intersect(cc$name, th_rows)
  neither <- setdiff(all_th, c(cc$name, th_rows))
  if (length(both)) cat("  BAD both:", both, "\n")
  if (length(neither)) cat("  BAD neither:", neither, "\n")
  for (i in seq_len(nrow(cc))) {
    r <- cc[i, ]
    sp <- tryCatch(set_prior(dens[[r$class]], class = r$class, group = r$group,
                             coef = r$coef, resp = r$resp, dpar = r$dpar,
                             nlpar = r$nlpar), error = function(e) e)
    if (inherits(sp, "error")) {
      cat("  BAD set_prior", r$name, r$class, ":", conditionMessage(sp), "\n")
      next
    }
    rr <- tryCatch(reached(des, sp), error = function(e) e)
    ok_reach <- !inherits(rr, "error") && r$name %in% rr
    # the slot must appear in the table too
    in_tab <- any(tab$class == r$class & tab$group == r$group &
                    tab$coef == r$coef & tab$resp == r$resp &
                    tab$dpar == r$dpar & tab$nlpar == r$nlpar)
    m <- msg_of(frmtmb:::resolve_priorlist(des, set_prior(
      "normal(0, 1)", class = "theta", coef = r$name)))
    ok_ref <- grepl(paste0("class = \"", r$class, "\""), m, fixed = TRUE)
    cat(sprintf("  %-9s %-6s grp=%-10s coef=%-14s resp=%-3s dpar=%-6s nlpar=%-3s reach=%s intable=%s refused=%s%s\n",
                r$name, r$class, r$group, r$coef, r$resp, r$dpar, r$nlpar,
                ok_reach, in_tab, ok_ref,
                if (inherits(rr, "error")) paste0(" [", conditionMessage(rr), "]")
                else if (!ok_reach) paste0(" [reached: ", paste(rr, collapse = ","), "]")
                else ""))
  }
  for (nm in th_rows) {
    m <- msg_of(frmtmb:::resolve_priorlist(des, set_prior(
      "normal(0, 1)", class = "theta", coef = nm)))
    cat(sprintf("  theta row %-10s resolves=%s %s\n", nm, !nzchar(m),
                substr(m, 1, 100)))
  }
  # the class-wide theta row resolves where it is listed, and is refused
  # where it is not
  m <- msg_of(frmtmb:::resolve_priorlist(des, set_prior("normal(0, 1)",
                                                        class = "theta")))
  cat("  classwide theta: ", if (nzchar(m)) substr(m, 1, 90) else "resolves",
      "\n")
  invisible()
}

probe("(1 | g1:g2)", bf(y ~ x + (1 | g2:g3)))
probe("(1 | g/g2)", bf(y ~ x + (1 | g/g2)))
probe("(1 + x || g)", bf(y ~ x + (1 + x || g)))
probe("gr(g, by = h)", bf(y ~ x + (1 + x | gr(g, by = h))))
probe("mm(g2, g3)", bf(y ~ x + (1 | mm(g2, g3))))
probe("gr(g, cov = A)", bf(y ~ x + (1 + x | gr(g, cov = A))), data2 = d2)
probe("gr(g, dist = student)", bf(y ~ x + (1 + x | gr(g, dist = "student"))))
probe("diag(0 + tf | g)", bf(y ~ x + diag(0 + tf | g)))
probe("homdiag(0 + tf | g)", bf(y ~ x + homdiag(0 + tf | g)))
probe("cs(0 + tf | g)", bf(y ~ x + cs(0 + tf | g)))
probe("homcs(0 + tf | g)", bf(y ~ x + homcs(0 + tf | g)))
probe("ar1(0 + tf | g)", bf(y ~ x + ar1(0 + tf | g)))
probe("hetar1(0 + tf | g)", bf(y ~ x + hetar1(0 + tf | g)))
probe("toep(0 + tf | g)", bf(y ~ x + toep(0 + tf | g)))
probe("ou(tn + 0 | g)", bf(y ~ x + ou(tn + 0 | g)))
probe("rr(0 + tf | g, d = 1)", bf(y ~ x + rr(0 + tf | g, d = 1)))
probe("gp(x) + gp(z)", bf(y ~ gp(x) + gp(z)))
probe("gp(x, by = h)", bf(y ~ gp(x, by = h)))
probe("gp(x, k = 6) hsgp", bf(y ~ gp(x, k = 6)))
probe("gp(x, z)", bf(y ~ gp(x, z)))
probe("s(x) + t2(x, z)", bf(y ~ s(x) + t2(x, z)))
probe("sigma ~ (1 | g)", bf(y ~ x + (1 | g), sigma ~ (1 | g)))
probe("sigma ~ (1 + x | g)", bf(y ~ x, sigma ~ x + (1 + x | g)))
probe("mv (1 | p | g)", bf(y ~ x + (1 | p | g)) + bf(y2 ~ x + (1 | p | g)) +
        set_rescor(FALSE))
probe("mv cs + rescor", bf(y ~ x + cs(0 + tf | g)) + bf(y2 ~ x) +
        bf(y3 ~ x) + set_rescor(TRUE))
probe("mv us per resp", bf(y ~ x + (1 + x | g)) + bf(y2 ~ x + (1 + x | g)) +
        set_rescor(TRUE))
probe("nonlinear", bf(yn ~ a * exp(k * x), a ~ 1 + (1 | g), k ~ 1 + (1 | g),
                      nl = TRUE))
probe("nonlinear shared ID", bf(yn ~ a * exp(k * x), a ~ 1 + (1 | q | g),
                                k ~ 1 + (1 | q | g), nl = TRUE))
probe("mixture mu1 (1 | g)", bf(y ~ 1, mu1 ~ x + (1 | g), mu2 ~ (1 | g)),
      family = mixture(gaussian(), gaussian()))
probe("categorical (1 | g)", bf(cat ~ x + (1 | g)), family = categorical())
probe("cumulative (1 + x | g)", bf(ord ~ x + (1 + x | g)),
      family = cumulative())
