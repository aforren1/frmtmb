# per fit: lane exact at its estimate against glmer's exact at glmer's
# estimate (short = glmer - lane, > 0 means the lane stopped below) and
# against base's exact; files: glmer rows, base rows (optional), lane rows
BEGIN { FS = "\t" }
$1 == "ROW" && $5 != "ERROR" {
  k = $2 " " $3; arm = $4
  for (i = 5; i < NF; i += 2) v[arm, k, $i] = $(i + 1)
  if (arm == "lane") keys[k] = 1
}
END {
  for (k in keys) {
    n++
    sh = (("glmer", k, "ref") in v) ? v["glmer", k, "ref"] - v["lane", k, "ref"] : "NA"
    sb = (("base", k, "ref") in v) ? v["base", k, "ref"] - v["lane", k, "ref"] : "NA"
    e = v["lane", k, "dref"] + 0; if (e < 0) e = -e; if (e > me) me = e
    if (sh != "NA" && sh > ms) { ms = sh; mk = k }
    if (sh != "NA" && sh > 1e-5) s5++
    if (sh != "NA" && sh > 1e-6) s6++
    if (v["lane", k, "warn"] != "-") { w++; print "  WARN", k, v["lane", k, "warn"] }
    if (v["lane", k, "rounds"] > 0) r++
    if (sh != "NA" && sh > 1e-6) printf "  SHORT %s glmer-lane %.3e base-lane %s rounds %s sgain %s sgrad %s dref %s\n", k, sh, sb, v["lane", k, "rounds"], v["lane", k, "sgain"], v["lane", k, "sgrad"], v["lane", k, "dref"]
    t += v["lane", k, "secs"]
  }
  printf "fits %d warned %d rounds>0 %d short>1e-6 %d short>1e-5 %d max short %.3e (%s) max|err| %.3e secs %.1f\n", n, w, r, s6, s5, ms, mk, me, t
}
