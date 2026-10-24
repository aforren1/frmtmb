#!/usr/bin/env bash
# Release review of 0.74.0: render the changed and new vignettes, 4 at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-release
export NOT_CRAN=true PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/review-rel074-log/render; mkdir -p $L
CO="C:/Users/adf44/source/r/rev-rel074-lib/core"
jobs_list=(
"vignettes/curve-inference.Rmd|"
"vignettes/frmtmb.Rmd|"
"vignettes/brms-migration.Rmd|"
"vignettes/inputs.Rmd|"
"vignettes/case-studies.Rmd|"
"extensions/frmtmb.latent/vignettes/latent.Rmd|"
"extensions/frmtmb.spatial/vignettes/spatial.Rmd|"
"extensions/frmtmb.survival/vignettes/royston-parmar.Rmd|"
"extensions/frmtmb.survival/vignettes/articles/spline-history.Rmd|"
"vignettes/frmtmb.Rmd|$CO"
"vignettes/curve-inference.Rmd|$CO"
)
i=0
for j in "${jobs_list[@]}"; do
  f="${j%%|*}"; lib="${j#*|}"
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 2; done
  i=$((i+1))
  if [ -n "$lib" ]; then
    ( "$RS" dev/review-rel074-render.R "$f" "$lib;C:/Users/adf44/AppData/Local/R/win-library/4.6" > $L/job$i-corelib.txt 2>&1 ) &
  else
    ( "$RS" dev/review-rel074-render.R "$f" > $L/job$i.txt 2>&1 ) &
  fi
done
wait
grep -h "^RENDER" $L/job*.txt
echo DONE
