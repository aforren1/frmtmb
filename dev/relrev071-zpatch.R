# brms_block_Z() of test-double-bar.R reads Z_<id>_<k>; a dpar's term is
# Z_<id>_<dpar>_<k> in brms's standata. Same layout, with the dpar infix.
brms_block_Z <- function(sd, id) {
  n <- sd$N
  L <- sd[[paste0("N_", id)]]
  nm <- grep(paste0("^Z_", id, "_"), names(sd), value = TRUE)
  k <- as.integer(regmatches(nm, regexpr("[0-9]+$", nm)))
  nm <- nm[order(k)]
  do.call(rbind, lapply(nm, function(z) {
    M <- matrix(0, n, L)
    M[cbind(seq_len(n), sd[[paste0("J_", id)]])] <- sd[[z]]
    M
  }))
}
