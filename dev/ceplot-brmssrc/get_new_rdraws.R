get_new_rdraws <- 
function (reframe, gf, rdraws, used_levels, old_levels, sample_new_levels, 
    draws = NULL) 
{
    snl_options <- c("uncertainty", "gaussian", "old_levels")
    sample_new_levels <- match.arg(sample_new_levels, snl_options)
    g <- unique(reframe$group)
    stopifnot(length(g) == 1)
    stopifnot(is.list(gf))
    used_by_per_level <- attr(used_levels, "by")
    old_by_per_level <- attr(old_levels, "by")
    new_levels <- setdiff(used_levels, old_levels)
    nranef <- nrow(reframe)
    nlevels <- length(old_levels)
    max_level <- nlevels
    out <- vector("list", length(gf))
    for (i in seq_along(gf)) {
        has_new_levels <- any(gf[[i]] > nlevels)
        if (has_new_levels) {
            new_indices <- sort(setdiff(gf[[i]], seq_len(nlevels)))
            out[[i]] <- matrix(NA, nrow(rdraws), nranef * length(new_indices))
            if (sample_new_levels == "uncertainty") {
                for (j in seq_along(new_indices)) {
                  if (length(old_by_per_level)) {
                    new_by <- used_by_per_level[used_levels == 
                      new_levels[j]]
                    possible_levels <- old_levels[old_by_per_level == 
                      new_by]
                    possible_levels <- which(old_levels %in% 
                      possible_levels)
                    sel_levels <- sample(possible_levels, NROW(rdraws), 
                      TRUE)
                  }
                  else {
                    sel_levels <- sample(seq_len(nlevels), NROW(rdraws), 
                      TRUE)
                  }
                  for (k in seq_len(nranef)) {
                    for (s in seq_rows(rdraws)) {
                      sel <- (sel_levels[s] - 1) * nranef + k
                      out[[i]][s, (j - 1) * nranef + k] <- rdraws[s, 
                        sel]
                    }
                  }
                }
            }
            else if (sample_new_levels == "old_levels") {
                for (j in seq_along(new_indices)) {
                  if (length(old_by_per_level)) {
                    new_by <- used_by_per_level[used_levels == 
                      new_levels[j]]
                    possible_levels <- old_levels[old_by_per_level == 
                      new_by]
                    possible_levels <- which(old_levels %in% 
                      possible_levels)
                    sel_level <- sample(possible_levels, 1)
                  }
                  else {
                    sel_level <- sample(seq_len(nlevels), 1)
                  }
                  for (k in seq_len(nranef)) {
                    sel <- (sel_level - 1) * nranef + k
                    out[[i]][, (j - 1) * nranef + k] <- rdraws[, 
                      sel]
                  }
                }
            }
            else if (sample_new_levels == "gaussian") {
                if (any(!reframe$dist %in% "gaussian")) {
                  stop2("Option sample_new_levels = 'gaussian' is not ", 
                    "available for non-gaussian group-level effects.")
                }
                for (j in seq_along(new_indices)) {
                  if (length(old_by_per_level)) {
                    new_by <- used_by_per_level[used_levels == 
                      new_levels[j]]
                    rnames <- as.vector(get_rnames(reframe, bylevels = new_by))
                  }
                  else {
                    rnames <- get_rnames(reframe)
                  }
                  sd_pars <- paste0("sd_", g, "__", rnames)
                  sd_draws <- prepare_draws(draws, sd_pars)
                  cor_type <- paste0("cor_", g)
                  cor_pars <- get_cornames(rnames, cor_type, 
                    brackets = FALSE)
                  cor_draws <- matrix(0, nrow(sd_draws), length(cor_pars))
                  for (k in seq_along(cor_pars)) {
                    if (cor_pars[k] %in% colnames(draws)) {
                      cor_draws[, k] <- prepare_draws(draws, 
                        cor_pars[k])
                    }
                  }
                  cov_matrix <- get_cov_matrix(sd_draws, cor_draws)
                  indices <- ((j - 1) * nranef + 1):(j * nranef)
                  out[[i]][, indices] <- t(apply(cov_matrix, 
                    1, rmulti_normal, n = 1, mu = rep(0, length(sd_pars))))
                }
            }
            max_level <- max_level + length(new_indices)
        }
        else {
            out[[i]] <- matrix(nrow = nrow(rdraws), ncol = 0)
        }
    }
    out <- do_call(cbind, out)
    structure(out, gf = gf, max_level = max_level)
}
