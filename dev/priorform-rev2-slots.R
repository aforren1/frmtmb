# Reviewer recheck round 1, lane wt-priorform: the duplicate-prior slot
# key and specificity precedence, brms 2.23.0 against the lane.
#   Rscript dev/priorform-rev2-slots.R brms|lane     seed 20260916
# Writes dev/priorform-rev2-slots-<mode>.txt
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (mode == "lane") "C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (mode == "brms") suppressMessages(library(brms)) else suppressMessages(library(frmtmb))
set.seed(20260916)
ng <- 10; n <- 200
d <- data.frame(g = factor(rep(1:ng, 20)), x = rnorm(n), z = rnorm(n),
                sigma = rnorm(n), a_x = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(ng)[d$g] + rnorm(n)
d$y1 <- d$y + rnorm(n); d$y2 <- rnorm(n) + d$z
d$yb <- pmin(pmax(plogis(d$x / 3 + rnorm(n, 0, .5)), .01), .99)
d$ynl <- 2 * exp(-0.3 * d$x) + rnorm(n, 0, .3)
out <- character(0)
say <- function(...) out <<- c(out, paste0(...))
P <- function(...) if (mode == "brms") brms::set_prior(...) else frmtmb::set_prior(...)
fam_g <- function() gaussian()
fam_b <- function() if (mode == "brms") brms::Beta() else frmtmb::Beta()
F <- function(txt) eval(parse(text = txt))
vp <- function(label, prior, form, fam = NULL, show = NULL) {
  r <- tryCatch({
    t <- if (is.null(fam)) validate_prior(prior, form, data = d) else
      validate_prior(prior, form, data = d, family = fam)
    t <- as.data.frame(t)
    if (mode == "brms") {
      # brms leaves an inheriting row empty; fill it from its parent
      # (coef row from the group row, group row from the class row) in
      # the same class, resp, dpar and nlpar, as brms applies it
      for (i in seq_len(nrow(t))) {
        if (nzchar(t$prior[i])) next
        same <- t$class == t$class[i] & t$resp == t$resp[i] &
          t$dpar == t$dpar[i] & t$nlpar == t$nlpar[i] & nzchar(t$prior)
        up <- if (nzchar(t$coef[i]) && nzchar(t$group[i]))
          which(same & t$group == t$group[i] & !nzchar(t$coef)) else integer(0)
        if (!length(up)) up <- which(same & !nzchar(t$group) & !nzchar(t$coef))
        if (length(up) && (nzchar(t$coef[i]) || nzchar(t$group[i])))
          t$prior[i] <- t$prior[up[1]]
      }
    }
    t$prior[!nzchar(t$prior)] <- "(flat)"
    if (!is.null(show)) {
      keep <- eval(show, t)
      sel <- t[keep, c("prior", "class", "coef", "group", "resp", "dpar", "nlpar")]
      sel <- sel[do.call(order, unname(as.list(sel[, -1]))), ]
      paste0("accept | ", paste(apply(sel, 1, function(r) paste(r[-1][nzchar(r[-1])], collapse = "/") |> paste0("=", r[1])), collapse = "; "))
    } else "accept"
  }, error = function(e) paste("REFUSE:", substr(gsub("\n", " ", conditionMessage(e)), 1, 110)))
  say(sprintf("%-34s %s", label, r))
}
# ---- slot keys: things brms accepts --------------------------------------
say("== slot keys")
f1 <- if (mode == "brms") bf(y ~ x + sigma + (1 | g), sigma ~ z) else bf(y ~ x + sigma + (1 | g), sigma ~ z)
vp("b dpar=sigma + b coef=sigma", P("normal(0,1)", class = "b", dpar = "sigma") +
     P("normal(0,2)", class = "b", coef = "sigma"), f1, fam_g())
vp("b coef=a_x + b nlpar... n/a", P("normal(0,1)", class = "b", coef = "a_x") +
     P("normal(0,2)", class = "b", coef = "x"), y ~ a_x + x, fam_g())
fmv <- if (mode == "brms") bf(y1 ~ x + (1 | g)) + bf(y2 ~ y1 + (1 | g)) + set_rescor(FALSE) else
  mvbf(bf(y1 ~ x + (1 | g)), bf(y2 ~ y1 + (1 | g)))
vp("mv b resp=y1 + b coef=y1 resp=y2", P("normal(0,1)", class = "b", resp = "y1") +
     P("normal(0,2)", class = "b", coef = "y1", resp = "y2"), fmv, if (mode == "brms") NULL else NULL)
vp("mv b resp=y2 + b coef=y1 resp=y2", P("normal(0,1)", class = "b", resp = "y2") +
     P("normal(0,2)", class = "b", coef = "y1", resp = "y2"), fmv)
vp("sd + sd group=g", P("normal(0,1)", class = "sd") + P("normal(0,2)", class = "sd", group = "g"),
   y ~ x + (1 | g), fam_g())
vp("sd + sd (same) refused?", P("normal(0,1)", class = "sd") + P("normal(0,2)", class = "sd"),
   y ~ x + (1 | g), fam_g())
vp("b coef='' + b (same)", P("normal(0,1)", class = "b", coef = "") + P("normal(0,2)", class = "b"),
   y ~ x, fam_g())
vp("Intercept + Intercept dpar=sigma", P("normal(0,1)", class = "Intercept") +
     P("normal(0,2)", class = "Intercept", dpar = "sigma"), bf(y ~ x, sigma ~ z), fam_g())
vp("nl b nlpar=a + b nlpar=a coef=x", P("normal(2,1)", nlpar = "a") +
     P("normal(0,2)", nlpar = "a", coef = "x"), bf(ynl ~ a * exp(-b * x), a ~ 1 + x, b ~ 1, nl = TRUE), fam_g())
vp("nl b nlpar=a coef=Intercept x2", P("normal(2,1)", nlpar = "a", coef = "Intercept") +
     P("normal(0,2)", nlpar = "a", coef = "Intercept"), bf(ynl ~ a * exp(-b * x), a ~ 1 + x, b ~ 1, nl = TRUE), fam_g())
vp("density then bounds same slot", P("normal(0,1)", coef = "x") + P("", coef = "x", lb = 0),
   y ~ x, fam_g())
# ---- precedence ------------------------------------------------------------
say("== precedence (table rows after validate)")
bshow <- quote(class == "b")
vp("b coef x THEN class b", P("normal(3,1)", coef = "x") + P("normal(0,1)", class = "b"), y ~ x + z, fam_g(), bshow)
vp("b class THEN coef x", P("normal(0,1)", class = "b") + P("normal(3,1)", coef = "x"), y ~ x + z, fam_g(), bshow)
sdshow <- quote(class == "sd")
fb <- bf(yb ~ x + (1 | g), phi ~ (1 | g))
vp("sd group g THEN sd", P("normal(0,5)", class = "sd", group = "g") + P("normal(0,1)", class = "sd"), fb, fam_b(), sdshow)
vp("sd dpar phi THEN sd", P("normal(0,5)", class = "sd", dpar = "phi") + P("normal(0,1)", class = "sd"), fb, fam_b(), sdshow)
vp("sd group g dpar phi THEN sd group g", P("normal(0,5)", class = "sd", group = "g", dpar = "phi") +
     P("normal(0,1)", class = "sd", group = "g"), fb, fam_b(), sdshow)
fc <- y ~ x + (x | g)
corshow <- quote(class %in% c("cor", "L"))
vp("cor group g THEN cor", P("lkj(4)", class = "cor", group = "g") + P("lkj(2)", class = "cor"), fc, fam_g(), corshow)
vp("b dpar sigma THEN b", P("normal(0,5)", class = "b", dpar = "sigma") + P("normal(0,1)", class = "b"),
   bf(y ~ x, sigma ~ z), fam_g(), bshow)
vp("nl coef x THEN nlpar a", P("normal(0,5)", nlpar = "a", coef = "x") + P("normal(2,1)", nlpar = "a"),
   bf(ynl ~ a * exp(-b * x), a ~ 1 + x, b ~ 1, nl = TRUE), fam_g(), bshow)
vp("mv sd group g THEN sd resp y1", P("normal(0,5)", class = "sd", group = "g") + P("normal(0,1)", class = "sd", resp = "y1"),
   fmv, NULL, sdshow)
vp("mv sd resp y1 THEN sd group g", P("normal(0,1)", class = "sd", resp = "y1") + P("normal(0,5)", class = "sd", group = "g"),
   fmv, NULL, sdshow)
vp("mv b coef y1 THEN b resp y2", P("normal(0,5)", class = "b", coef = "y1") + P("normal(0,1)", class = "b", resp = "y2"),
   fmv, NULL, bshow)
vp("Intercept dpar sigma THEN Intercept", P("normal(0,5)", class = "Intercept", dpar = "sigma") +
     P("normal(0,1)", class = "Intercept"), bf(y ~ x, sigma ~ z), fam_g(), quote(class == "Intercept"))
writeLines(out, sprintf("C:/Users/adf44/source/r/frmtmb-wt-priorform/dev/priorform-rev2-slots-%s.txt", mode))
