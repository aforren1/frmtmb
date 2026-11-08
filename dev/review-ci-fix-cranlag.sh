#!/usr/bin/env bash
# Short check of the fix round, item 1: run the revised
# dev/release/cran-lag.R on the cases that failed open, and on new ones.
# Scratch roots (copies of the DESCRIPTION files and the script) live in
# C:/Users/adf44/source/r/rev-ci-src; nothing in the worktree changes.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-ci
SRC=/c/Users/adf44/source/r/rev-ci-src
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
OUT=$WT/dev/review-ci-fix-cranlag.log
USER_LIB="C:/Users/adf44/AppData/Local/R/win-library/4.6"
STACK="C:/Users/adf44/source/r/rellib-r18;$USER_LIB"
: > "$OUT"

# a copy of the repository's DESCRIPTION files and the script
mkroot() {
  local r=$SRC/$1
  rm -rf "$r"; mkdir -p "$r/dev/release"
  cp "$WT/DESCRIPTION" "$r/"
  for d in "$WT"/extensions/*/; do
    p=$(basename "$d"); mkdir -p "$r/extensions/$p"
    cp "$d/DESCRIPTION" "$r/extensions/$p/"
  done
  cp "$WT/dev/release/cran-lag.R" "$r/dev/release/"
  echo "$r"
}
run() {
  local name=$1; shift
  echo "=== $name: cran-lag.R $*" >> "$OUT"
  "$@" >> "$OUT.tmp" 2>&1
  local st=$?
  grep -v -E '^ +[A-Za-z0-9.]+ +current' "$OUT.tmp" | tail -40 >> "$OUT"
  rm -f "$OUT.tmp"
  echo "--- $name exit $st" >> "$OUT"
}

S=$WT/dev/release/cran-lag.R
run A_bad_mirror "$R" "$S" "$USER_LIB" "http://127.0.0.1:9"
run B_offline env HTTP_PROXY=http://127.0.0.1:9 \
  HTTPS_PROXY=http://127.0.0.1:9 http_proxy=http://127.0.0.1:9 \
  https_proxy=http://127.0.0.1:9 "$R" "$S" "$USER_LIB"

# C: a Suggests package that no repository carries and nothing installs
r=$(mkroot root-ghost)
sed -i 's/^    withr$/    withr,\n    ghostpkgreviewci/' "$r/DESCRIPTION"
grep -n ghostpkgreviewci "$r/DESCRIPTION" >> "$OUT"
run C_ghost_suggests "$R" "$r/dev/release/cran-lag.R" "$STACK"

# D: the stack with rellib-r18 (drmTMB lives there)
run D_stack "$R" "$S" "$STACK"

# E: CRAN up, the r-universe repository down
r=$(mkroot root-universe)
for f in "$r/DESCRIPTION" "$r"/extensions/*/DESCRIPTION; do
  sed -i 's#https://kaskr.r-universe.dev#http://127.0.0.1:9#' "$f"
done
grep -rn "127.0.0.1" "$r" --include=DESCRIPTION | head -3 >> "$OUT"
run E_universe_down "$R" "$r/dev/release/cran-lag.R" "$STACK"

# F: malformed library arguments
run F1_missing_dir "$R" "$S" "C:/no/such/lib;$USER_LIB"
run F2_double_sep "$R" "$S" "$USER_LIB;;C:/Users/adf44/source/r/rellib-r18"
run F3_not_a_library "$R" "$S" "C:/Users/adf44/source/r/frmtmb-wt-ci/dev"
run F4_lone_sep "$R" "$S" ";"

# G: a file:// mirror that holds exactly the stack (all current), with
# no Windows binary index; then the same mirror cut to half its entries
"$R" "$WT/dev/ci-lagrepo.R" "$SRC/lagrepo-full" "$STACK" >> "$OUT" 2>&1
run G1_all_current "$R" "$S" "$STACK" "file:///C:/Users/adf44/source/r/rev-ci-src/lagrepo-full"
"$R" -e 'p <- read.dcf("C:/Users/adf44/source/r/rev-ci-src/lagrepo-full/src/contrib/PACKAGES"); dir.create("C:/Users/adf44/source/r/rev-ci-src/lagrepo-half/src/contrib", recursive = TRUE, showWarnings = FALSE); keep <- seq_len(nrow(p)) %% 2 == 0; write.dcf(p[keep, , drop = FALSE], "C:/Users/adf44/source/r/rev-ci-src/lagrepo-half/src/contrib/PACKAGES"); cat("half:", sum(keep), "of", nrow(p), "\n")' >> "$OUT" 2>&1
run G2_partial_index "$R" "$S" "$STACK" "file:///C:/Users/adf44/source/r/rev-ci-src/lagrepo-half"
echo done >> "$OUT"
