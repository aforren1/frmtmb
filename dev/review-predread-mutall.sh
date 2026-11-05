#!/usr/bin/env bash
# All mutants, one after another, each with the audit shapes it reaches.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
run() { bash dev/review-predread-mutrun.sh "$1" "$2" >> dev/review-predread-out/mutall.log 2>&1; }
run cs_se cs,cs_cum
run cs_se_rev cs,cs_cum
run cs_factor_nd cs_cum
run me_var_double me
run me_var_drop me
run me_half_sd me
run me_gr_perrow me
run link_identity cumulative,cs
run emm_me_design me
run me_store_off me
run curve_check_ratio cumulative
run resp_read_off plain
run emm_cs_design cs,cs_cum
run me_share me
