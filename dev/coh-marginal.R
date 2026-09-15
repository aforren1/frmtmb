## Round 1 punch: is the low rungs' offset a bias or a different
## estimand?
##
## The rungs are PAIRED, since every rung of a replicate saw the same
## data at the same n, so the difference of two rungs is measured far
## more sharply than either against 0.5. And a logit contrast is not
## collapsible: a model that leaves variance in the linear predictor
## unmodeled is consistent for the POPULATION-AVERAGED contrast, which
## is smaller than the conditional one by about
##
##     1 / sqrt(1 + 0.346 V),  V the omitted variance
##
## This computes both sides: the paired differences from
## dev/coh-recovery-main.tsv, and the prediction from the generator's
## own constants in dev/coh-sim.R with nothing fitted.
##
## Run: Rscript dev/coh-marginal.R
src <- readLines("dev/coh-summarize.R")
stop_at <- grep("^truth <- ", src)[1L]
eval(parse(text = paste(src[seq_len(stop_at - 1L)], collapse = "\n")))
source("dev/coh-sim.R")

d <- read_tsv_kv("dev/coh-recovery-main.tsv")
keep <- names(which(table(d$seed) == length(unique(d$rung))))
d <- d[as.character(d$seed) %in% keep, ]
w <- reshape(d[, c("seed", "rung", "est")], idvar = "seed",
             timevar = "rung", direction = "wide")
names(w) <- sub("^est[.]", "", names(w))
cat("complete replicates:", nrow(w), "\n\n")

pair <- function(a, b, lab) {
  x <- w[[a]] - w[[b]]
  se <- stats::sd(x) / sqrt(length(x))
  cat(sprintf("%-28s %9.5f  se %7.5f  t %7.2f\n", lab, mean(x), se,
              mean(x) / se))
  invisible(mean(x))
}
cat("against the correct model:\n")
b_cond <- pair("cond", "full", "cond minus full")
b_smooth <- pair("smooth", "full", "smooth minus full")
b_id <- pair("id", "full", "id minus full")
pair("idcond", "full", "idcond minus full")
cat("\nrung to rung, which term moves it:\n")
pair("smooth", "id", "dropping (1 | id)")
pair("cond", "smooth", "dropping the smooth")
pair("id", "full", "dropping (1 | id:cond)")

## The omitted variance, from the truth and not from any fit. The
## design is balanced over the same frequency grid in every cell, so
## the bump contributes its own variance over that grid.
tr <- coupling_truth
fr <- seq_len(tr$n_freq) / tr$n_freq
bump <- coupling_fbump(fr)
v_bump <- mean((bump - mean(bump))^2)
v_id <- tr$sd_id^2
v_ic <- tr$sd_idcond^2
att <- function(V) 1 / sqrt(1 + 0.346 * V)
cat(sprintf("\nbump variance over the design grid: %.5f\n", v_bump))
cat(sprintf("V(cond) = %.4f   V(smooth) = %.4f   V(id) = %.4f\n",
            v_bump + v_id + v_ic, v_id + v_ic, v_ic))
tab <- rbind(
  c(v_bump + v_id + v_ic, tr$b_cond * (att(v_bump + v_id + v_ic) - 1),
    b_cond),
  c(v_id + v_ic, tr$b_cond * (att(v_id + v_ic) - 1), b_smooth),
  c(v_ic, tr$b_cond * (att(v_ic) - 1), b_id))
dimnames(tab) <- list(c("cond", "smooth", "id"),
                      c("V", "predicted", "observed"))
print(round(tab, 5))
cat("\nresidual, observed minus predicted:\n")
print(round(tab[, "observed"] - tab[, "predicted"], 5))

## The rung-to-rung arithmetic, to check that the one residual the
## model misses is carried by the `id` rung and not spread over the
## ladder: what the attenuation predicts for dropping (1 | id) is the
## difference of two predictions.
p_smooth <- tr$b_cond * (att(v_id + v_ic) - 1)
p_id <- tr$b_cond * (att(v_ic) - 1)
cat(sprintf(paste0("\ndropping (1 | id): predicted %.5f, observed",
                   " %.5f, difference %.5f\n"),
            p_smooth - p_id, mean(w$smooth - w$id),
            mean(w$smooth - w$id) - (p_smooth - p_id)))
cat(sprintf("the id rung's own residual is %.5f\n",
            b_id - p_id))

## And that last line is an IDENTITY, not a second check: the paired
## means and the predictions both subtract across rungs, so residuals
## subtract exactly. Printed because "the residual sits on one rung"
## means something only once the additivity is known to be arithmetic.
r_diff <- mean(w$smooth - w$id) - (p_smooth - p_id)
r_smooth <- b_smooth - p_smooth
r_id <- b_id - p_id
cat(sprintf(paste0("r(smooth - id) + r(id) - r(smooth) = %.3g;",
                   "  all.equal: %s\n"),
            r_diff + r_id - r_smooth,
            isTRUE(all.equal(r_diff, r_smooth - r_id))))

## The same arithmetic in the NULL arm, where the truth has no
## id-by-condition component, so `id` is the correct model and the
## omitted variance of the two lowest rungs is smaller by that amount.
dn <- read_tsv_kv("dev/coh-recovery-null.tsv")
kn <- names(which(table(dn$seed) == length(unique(dn$rung))))
dn <- dn[as.character(dn$seed) %in% kn, ]
wn <- reshape(dn[, c("seed", "rung", "est")], idvar = "seed",
              timevar = "rung", direction = "wide")
names(wn) <- sub("^est[.]", "", names(wn))
cat(sprintf("\nnull arm, %d replicates, paired against `id`:\n",
            nrow(wn)))
for (rg in c("cond", "smooth")) {
  V <- if (identical(rg, "cond")) v_bump + v_id else v_id
  x <- wn[[rg]] - wn$id
  se <- stats::sd(x) / sqrt(length(x))
  cat(sprintf("%-8s V %.5f  predicted %.5f  observed %.5f  se %.5f\n",
              rg, V, tr$b_cond * (att(V) - att(0)), mean(x), se))
}
