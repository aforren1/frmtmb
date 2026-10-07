#!/usr/bin/env bash
# Build a copy of R 4.6.1 whose BLAS (and optionally LAPACK) is a given
# OpenBLAS release, to see on this Windows machine the rounding of the
# Ubuntu runners. Generalizes dev/cifix-openblas.sh, which was fixed to
# 0.3.26, BLAS only, and to a worktree that no longer exists.
#
#   bash dev/ciharden-openblas.sh <version> <blas|lapack>
#
#   0.3.26 is what ubuntu-24.04 installs (libopenblas0-pthread
#   0.3.26+ds-1ubuntu0.1); 0.3.32 is what ubuntu-26.04 installs
#   (0.3.32+ds-5), checked on packages.ubuntu.com on 2026-10-06.
#
# "blas" keeps R's reference LAPACK, which calls OpenBLAS for its BLAS.
# "lapack" also routes every LAPACK routine that OpenBLAS exports to
# OpenBLAS, as on the runner, where the openblas-pthread alternative
# provides both libblas.so.3 and liblapack.so.3. OpenBLAS replaces some
# LAPACK routines (getrf, potrf, trtri, lauum, ...) with its own
# blocked, threaded code, so this arm rounds differently from "blas".
#
# R's DLLs import by name, so a forwarding DLL that exports every
# symbol of the original and points each one at libopenblas.dll, or
# at the renamed original when OpenBLAS lacks it (dgemmtr_ in older
# releases, R's Fortran module symbols), is all it takes. The installed
# R is not touched; the copy lands in dev/ciharden-out/Rob<ver>-<mode>.
#
# Threads: OpenBLAS reads OPENBLAS_NUM_THREADS (or OMP_NUM_THREADS)
# at load. The runners have 4 cores and set neither, so 4 is the
# runner's configuration.
set -eu
VER="$1"; MODE="$2"
case "$MODE" in blas|lapack) ;; *) echo "mode: blas or lapack"; exit 1;; esac
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/dev/ciharden-out"
DEST="$OUT/Rob$VER-$MODE"
export PATH="/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
mkdir -p "$OUT/openblas-$VER" "$OUT/shim-$VER-$MODE"
cd "$OUT/openblas-$VER"
if [ ! -f bin/libopenblas.dll ]; then
  gh release download "v$VER" -R OpenMathLib/OpenBLAS \
    -p "OpenBLAS-$VER-x64.zip" --clobber
  unzip -o -q "OpenBLAS-$VER-x64.zip"
fi
[ -f bin/libopenblas.dll ] || { echo "no bin/libopenblas.dll"; exit 1; }
SH="$OUT/shim-$VER-$MODE"
exports() {
  objdump -p "$1" |
    sed -n '/\[Ordinal\/Name Pointer\] Table/,/^$/p' |
    awk '{print $NF}' | tr -d '\r' | grep -E '^[A-Za-z_][A-Za-z0-9_]*$' |
    grep -v '^Name$' | sort -u
}
exports bin/libopenblas.dll > "$SH/openblas.txt"
rm -rf "$DEST"
cp -r "/c/Program Files/R/R-4.6.1" "$DEST"
X="$DEST/bin/x64"
cp bin/libopenblas.dll "$X/libopenblas.dll"
# forward <dll stem> <renamed original stem>
forward() {
  local stem="$1" ref="$2"
  cp "$X/$stem.dll" "$X/$ref.dll"
  exports "$X/$ref.dll" > "$SH/$stem-exports.txt"
  {
    echo "LIBRARY $stem.dll"
    echo "EXPORTS"
    while read -r s; do
      if grep -qx "$s" "$SH/openblas.txt"; then
        echo "  $s = libopenblas.$s"
      else
        echo "  $s = $ref.$s"
      fi
    done < "$SH/$stem-exports.txt"
  } > "$SH/$stem.def"
  echo "void ${stem}_shim_dummy(void) {}" > "$SH/$stem.c"
  gcc -shared -o "$SH/$stem.dll" "$SH/$stem.c" "$SH/$stem.def"
  cp "$SH/$stem.dll" "$X/$stem.dll"
  echo "$stem: $(grep -c 'libopenblas[.]' "$SH/$stem.def") to OpenBLAS," \
    "$(grep -c "= $ref[.]" "$SH/$stem.def") to $ref"
}
forward Rblas Rblasref
[ "$MODE" = lapack ] && forward Rlapack Rlapackref
"$X/Rscript.exe" "$ROOT/dev/ciharden-blascheck.R"
