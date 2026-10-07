#!/usr/bin/env bash
# Build a copy of R 4.6.1 whose BLAS is OpenBLAS 0.3.26, the version the
# Ubuntu runners install, so that a Linux-only rounding failure can be
# seen on this Windows machine. Everything lands in dev/setier-out/,
# which dev/.gitignore keeps out of the repository. The installed R is
# not touched.
#
# R's Rblas.dll exports dgemmtr_ and zgemmtr_ (LAPACK 3.12), which
# OpenBLAS 0.3.26 lacks, so Rlapack.dll will not load against a plain
# copy of libopenblas.dll. A forwarding Rblas.dll sends every other
# symbol to libopenblas.dll and those two to the reference BLAS, kept
# as Rblasref.dll. LAPACK stays R's reference build, calling OpenBLAS
# for its BLAS.
#
# Check: dev/cifix-blascheck.R times a 2000 x 2000 product at about
# 0.05 s against 4.4 s with the reference BLAS.
set -eu
ROOT=/c/Users/adf44/source/r/frmtmb-wt-setier
OUT="$ROOT/dev/setier-out"
export PATH="/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
mkdir -p "$OUT/openblas" "$OUT/shim"
cd "$OUT/openblas"
gh release download v0.3.26 -R OpenMathLib/OpenBLAS \
  -p "OpenBLAS-0.3.26-x64.zip" --clobber
unzip -o -q OpenBLAS-0.3.26-x64.zip
cd "$OUT"
rm -rf Rob
cp -r "/c/Program Files/R/R-4.6.1" Rob
X="$OUT/Rob/bin/x64"
cp "$X/Rblas.dll" "$X/Rblasref.dll"
cp openblas/bin/libopenblas.dll "$X/libopenblas.dll"
objdump -p "$X/Rblasref.dll" |
  sed -n '/\[Ordinal\/Name Pointer\] Table/,/^$/p' |
  awk '{print $NF}' | tr -d '\r' | grep -E '^[a-z][a-z0-9_]*$' |
  sort -u > shim/exports.txt
{
  echo "LIBRARY Rblas.dll"
  echo "EXPORTS"
  grep -v "gemmtr_" shim/exports.txt |
    awk '{print "  " $1 " = libopenblas." $1}'
  echo "  dgemmtr_ = Rblasref.dgemmtr_"
  echo "  zgemmtr_ = Rblasref.zgemmtr_"
} > shim/rblas.def
echo 'void rblas_shim_dummy(void) {}' > shim/shim.c
gcc -shared -o shim/Rblas.dll shim/shim.c shim/rblas.def
cp shim/Rblas.dll "$X/Rblas.dll"
"$X/Rscript.exe" "$(cygpath -w "$ROOT/dev/cifix-blascheck.R")"
