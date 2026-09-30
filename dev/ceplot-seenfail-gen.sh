# Lane ceplot: the failing expectations of each new test on the base
# build, from dev/ceplot-log/seen-failing-*.txt, one line each.
cd /c/Users/adf44/source/r/frmtmb-wt-ceplot
for f in dev/ceplot-log/seen-failing-*.txt; do
  echo "$(basename "$f"):"
  grep -a -A2 "^── [0-9]*\. \(Failure\|Error\)" "$f" | tr -d '\r' |
    grep -a -v "^--$\|^<frmtmb" | sed 's/^/  /' | cut -c1-180
done
