r <- readRDS("dev/review-equiv-scratch/results.rds")
sel <- commandArgs(trailingOnly = TRUE)[1]
a <- if (identical(sel, "all")) r else r[r$agree, ]
for (i in seq_len(nrow(a))) {
  cat(sprintf("%s | %s x %s | %s | %s | %d | %s | c=%.2g g=%.2g pts=%s | L2: %s | L4: %s %s | L1err: %s | frmerr: %s | brmserr: %s\n",
              a$id[i], a$feature_a[i], a$feature_b[i], a$status[i],
              a$level[i], a$draws[i], a$verdict[i], a$l1_const[i],
              a$l1_grad[i], a$l1_points[i], substr(a$l2_what[i], 1, 160),
              a$l4_route[i], a$l4_z[i], substr(a$l1_error[i], 1, 80),
              substr(a$frm_error[i], 1, 120), substr(a$brms_error[i], 1, 120)))
}
