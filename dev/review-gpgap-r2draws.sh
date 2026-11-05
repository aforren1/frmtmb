#!/usr/bin/env bash
# fix round 2, F1, F2, item 2: the installed build's drawn route at the settings of
# round 1's table, and brms at n 40.
set -u
R_PARALLEL_PORT=13011 bash dev/review-gpgap-run.sh dev/review-gpgap-draw.R r2-frm-11.txt frm 11 3000 12 &
R_PARALLEL_PORT=13057 bash dev/review-gpgap-run.sh dev/review-gpgap-draw.R r2-frm-22.txt frm 22 3000 12 &
R_PARALLEL_PORT=13103 bash dev/review-gpgap-run.sh dev/review-gpgap-draw.R r2-fit-11.txt frm_fit 11 3000 12 &
R_PARALLEL_PORT=13151 bash dev/review-gpgap-run.sh dev/review-gpgap-draw40.R r2-brms40-6601.txt brms 6601 3000 12 &
wait
