# summarize quadseam-review-bnd.R rows per design
/^ROW/ { des=$2; sing=$10; split($0, p, "[|]"); qb = (p[2] ~ /bmsg 1/); lb = (p[3] ~ /bmsg 1/); qn = (p[2] ~ /seNaN TRUE/)
  n[des]++; if (sing=="TRUE") gs[des]++; if (qb) q[des]++; if (lb) l[des]++
  if (qb && sing=="FALSE") { qf[des]++; print "  quad msg, glmer not singular:", $0 }
  if (lb && sing=="FALSE") lf[des]++
  if (qb != lb) { dis[des]++; print "  quad != lap:", $0 }
  if (sing=="TRUE" && !qb) { miss[des]++; print "  glmer singular, no quad msg:", $0 }
  if (qn && !qb) nanq[des]++ }
END { for (d in n) printf "%-8s n %d glmer_sing %d quad_msg %d lap_msg %d quad_msg_not_sing %d lap_msg_not_sing %d quad!=lap %d sing_no_quadmsg %d nanSE_no_msg %d\n", d, n[d], gs[d], q[d], l[d], qf[d], lf[d], dis[d], miss[d], nanq[d] }
