# join ROW lines of base, lane, glmer arms by design and seed
BEGIN { FS = "\t" }
$1 == "ROW" && $5 != "ERROR" {
  k = $2 " " $3; arm = $4
  for (i = 5; i < NF; i += 2) v[arm, k, $i] = $(i + 1)
  keys[k] = 1
}
$1 == "ROW" && $5 == "ERROR" { err[$4, $2 " " $3] = $6; keys[$2 " " $3] = 1 }
END {
  printf "%-12s %12s %12s %12s %12s %12s %5s %-40s %-40s\n", "design seed", "dref_base", "dref_lane", "exL-exB", "exG-exL", "sdL-sdB", "rnd", "warn_lane(base)", "msg_lane(base)"
  for (k in keys) {
    if ((("lane", k) in err) || (("base", k) in err)) { print k, "ERROR", err["base", k], err["lane", k]; continue }
    eg = (("glmer", k, "ref") in v) ? sprintf("%12.3e", v["glmer", k, "ref"] - v["lane", k, "ref"]) : sprintf("%12s", "-")
    printf "%-12s %12.3e %12.3e %12.3e %s %12.3e %5s %-40s %-40s\n", k, v["base", k, "dref"], v["lane", k, "dref"], v["lane", k, "ref"] - v["base", k, "ref"], eg, v["lane", k, "sd"] - v["base", k, "sd"], v["lane", k, "rounds"], v["lane", k, "warn"] "(" v["base", k, "warn"] ")", v["lane", k, "msg"] "(" v["base", k, "msg"] ")"
  }
}
