#!/usr/bin/env bash
# Merge the three lanes of the round opened 2026-10-28 (after 0.77.0)
# onto the release tree, file by file with git merge-file, in the order
# equiv, meme, rescorpp. equiv changes no tracked file (its survey is
# dev/ only); meme and rescorpp change different code files and meet in
# the two NEWS files only. Every lane and the release tree are on
# ac4906bd (frmtmb 0.77.0 plus one test fix, pushed), so each file
# starts from its ac4906bd blob and each lane's change is merged in
# three-way: the release file (ours), the ac4906bd blob (base) and the
# lane's file (theirs), LF-normalized. It touches no index and makes no
# git write: blobs are read with `git show`. Adapted from
# dev/rel077-merge.sh.
#
#   bash dev/rel078-merge.sh
#
# Rerunnable: every run first restores each file any lane changed to its
# HEAD blob and removes the untracked lane files an earlier run copied
# ($LOG/copied-<lane>.lst), then merges the lanes in order and copies
# their untracked files again. The hand resolutions are applied by
# dev/rel078-resolve.R (which sources dev/rel078-release-fix.R), the
# bin-2 tier files are regenerated from the adapter spec
# (dev/brmsport-bin2-gen.R), and the NEWS files are written by
# dev/rel078-news-merge.R, all run at the end. Last, every file the
# merge wrote is stored with CRLF line endings, the working copy's
# convention (core.autocrlf=true; the index holds LF either way).
#
# Generated files are left out of the textual merge, because their
# generators run after it: man/*.Rd and NAMESPACE (roxygen),
# frmtmb.sample's copy of helper-brms-suite.R, R/links-brms.R
# (dev/famlink-gen-links.R; the run checks that the regenerated file is
# the HEAD blob) and the bin-2 tier files under dev/brmsport-bin2/
# (dev/brmsport-bin2-gen.R). Untracked lane files are copied, except
# generated Rd files, Stan caches, the brms suite, OpenBLAS builds of R,
# the reviewers' output folders that dev/.gitignore lets through (a
# dev/review-*-out*/ folder: equiv's dev/review-equiv-out2/, its rerun
# results) and job lists under dev/release/. The reviewers' own scripts,
# the loose dev/review-meme-* and dev/review-rescorpp-* files and the 26
# scripts of equiv's dev/review-equiv-scratch/ (its logs, .rds and .txt
# files are gitignored), are copied, as the reviews themselves are: each
# review names them as its evidence (the 0.78.0 release review, finding
# 1, for the equiv ones).
set -u
# no git call here takes the index lock to refresh stat data, in this
# tree or in a lane's (the 0.77.0 release review, finding 4)
export GIT_OPTIONAL_LOCKS=0
BASE=ac4906bd
HEADC=ac4906bd
REL=${REL078_ROOT:-C:/Users/adf44/source/r/frmtmb-wt-rel078}
G=${REL078_GIT:-$REL}
export REL078_ROOT="$REL" REL078_GIT="$G"
LOG=$REL/dev/rel078-log
W=$LOG/merge-work
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LANES="equiv meme rescorpp"
export REL078_LANES="$LANES"
cd "$REL" || exit 1
[ "$(git -C "$G" rev-parse --short=8 HEAD)" = "$HEADC" ] ||
  { echo "release tree is not on $HEADC"; exit 1; }
mkdir -p "$LOG"; rm -rf "$W"; mkdir -p "$W/keep"
excl='^(man/.*\.Rd|NAMESPACE|R/links-brms\.R|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R|dev/brmsport-bin2/.*\.R)$'
uskip='^(man/.*\.Rd|extensions/[^/]*/man/.*\.Rd|.*stan-cache.*|dev/ciharden-out/.*|dev/brms-suite/.*|dev/review-[^/]*-out[^/]*/.*|dev/release/[^/]*\.jobs)$'
echo "lanes: $LANES (base $BASE, release HEAD $HEADC)"

: > "$W/all"
for l in $LANES; do
  src=/c/Users/adf44/source/r/frmtmb-wt-$l
  [ "$(git -C "$src" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "lane $l is not on $BASE"; exit 1; }
  git -C "$src" diff --name-status "$BASE" 2>/dev/null > "$W/$l.status"
  cut -f2 "$W/$l.status" >> "$W/all"
done
sort -u "$W/all" > "$W/all-files"
grep -vE "$excl" "$W/all-files" > "$W/files"
echo "files changed or deleted by some lane: $(wc -l < "$W/all-files")," \
     "without generated ones: $(wc -l < "$W/files")"

# restore every file any lane changed, generated ones too, so that a
# rerun starts from the release tree's own blob
while IFS= read -r f; do
  mkdir -p "$(dirname "$f")"
  git -C "$G" show "HEAD:$f" | tr -d '\r' > "$f"
done < "$W/all-files"
for l in $LANES; do
  m="$LOG/copied-$l.lst"
  [ -f "$m" ] || continue
  while IFS= read -r f; do
    git -C "$G" ls-files --error-unmatch "$f" > /dev/null 2>&1 && continue
    rm -f "$f"
  done < "$m"
  rm -f "$m"
done

: > "$W/deleted"
for l in $LANES; do
  src=/c/Users/adf44/source/r/frmtmb-wt-$l
  echo "== lane $l"
  nc=0; nk=0; nd=0
  while IFS=$'\t' read -r st f; do
    case "$st" in
      D)
        git -C "$G" show "$BASE:$f" | tr -d '\r' > "$W/base"
        if [ ! -e "$f" ]; then
          nd=$((nd + 1))
        elif tr -d '\r' < "$f" | cmp -s - "$W/base"; then
          rm -f "$f"; nd=$((nd + 1)); printf '%s\n' "$f" >> "$W/deleted"
        else
          echo "   CONFLICT deletion of a file an earlier lane changed: $f"
          nk=$((nk + 1))
        fi
        continue;;
    esac
    echo "$f" | grep -qE "$excl" && continue
    case "$st" in
      M)
        if [ ! -e "$f" ]; then
          echo "   CONFLICT change to a file an earlier lane deleted: $f"
          nk=$((nk + 1)); continue
        fi
        git -C "$G" show "$BASE:$f" | tr -d '\r' > "$W/base"
        tr -d '\r' < "$src/$f" > "$W/theirs"
        tr -d '\r' < "$f" > "$W/ours"
        k="$W/keep/$l.$(echo "$f" | tr '/' '_')"
        cp "$W/ours" "$k.ours"; cp "$W/base" "$k.base"
        cp "$W/theirs" "$k.theirs"
        if git merge-file -L release -L base -L "$l" \
             "$W/ours" "$W/base" "$W/theirs"; then
          rm -f "$k.ours" "$k.base" "$k.theirs"; nc=$((nc + 1))
        else
          echo "   CONFLICT($?) $f"; nk=$((nk + 1))
        fi
        cp "$W/ours" "$f";;
      *) echo "   UNHANDLED $st $f";;
    esac
  done < "$W/$l.status"
  echo "   merged clean: $nc, with conflicts: $nk, deleted: $nd"
  lst="$W/$l.untracked"; : > "$lst"; clash=0; skip=0
  while IFS= read -r f; do
    if echo "$f" | grep -qE "$uskip"; then skip=$((skip + 1)); continue; fi
    if [ -e "$f" ]; then
      echo "   CLASH untracked: $f"; clash=$((clash + 1)); continue
    fi
    printf '%s\n' "$f" >> "$lst"
  done < <(git -C "$src" ls-files --others --exclude-standard)
  ( cd "$src" && tar -cf - -T "$lst" ) | tar -xpf -
  # untracked copies keep the lane's line endings; store them as LF
  # here (the resolver reads LF), CRLF at the end
  while IFS= read -r f; do
    case "$f" in *.R|*.r|*.md|*.Rmd|*.sh|*.tsv|*.yaml|*.yml|*.txt|*.jobs|*.ps1|*.cmd|*.patch|*.sed|DESCRIPTION|NAMESPACE|.Rbuildignore)
      if grep -q $'\r' "$f"; then sed -i 's/\r$//' "$f"; fi;;
    esac
  done < "$lst"
  cp "$lst" "$LOG/copied-$l.lst"
  echo "   untracked copied: $(grep -c '' "$lst"), skipped: $skip," \
       "clashes: $clash"
done
echo "MERGE DONE (textual); conflict inputs kept in $W/keep"
grep -ln '^<<<<<<< ' $(grep -vxFf "$W/deleted" "$W/files" |
  while read -r f; do [ -e "$f" ] && echo "$f"; done) 2>/dev/null |
  sed 's/^/   conflicted after the textual merge: /'
[ "${REL078_TEXTUAL_ONLY:-}" = true ] && exit 0

# R/links-brms.R from the generator, which no lane changed; it reads brms
# from the user library, and must give the HEAD blob back
"$RS" dev/famlink-gen-links.R > "$LOG/links-gen.txt" 2>&1 || exit 1
if git -C "$G" show HEAD:R/links-brms.R | tr -d '\r' |
     cmp -s - <(tr -d '\r' < R/links-brms.R); then
  echo "R/links-brms.R regenerated: equal to the HEAD blob"
else
  echo "R/links-brms.R regenerated: DIFFERS from the HEAD blob"; exit 1
fi
"$RS" dev/rel078-resolve.R || exit 1
# the bin-2 tier files from the adapter spec, which no lane changed and
# release fix R2 edits (a pass note); generated, never hand-merged. The
# generator reads the verdicts, which are the HEAD blob's until the
# ledger is regenerated, so here the files come back as the HEAD blobs;
# the bin-2 pipeline (dev/brmsport-bin2-all.sh) then writes R2's note
# into the posterior_predict tier file.
"$RS" dev/brmsport-bin2-gen.R > "$LOG/bin2-gen.txt" 2>&1 ||
  { cat "$LOG/bin2-gen.txt"; exit 1; }
echo "bin-2 tier files regenerated ($(grep -c '' "$LOG/bin2-gen.txt") log lines)"
nb=0; nbd=0
for f in $(git -C "$G" ls-files 'dev/brmsport-bin2/*.R'); do
  nb=$((nb + 1))
  git -C "$G" show "HEAD:$f" | tr -d '\r' | cmp -s - <(tr -d '\r' < "$f") ||
    { nbd=$((nbd + 1)); echo "   bin-2 tier file differs from HEAD: $f"; }
done
echo "bin-2 tier files: $nb, differ from the HEAD blob: $nbd"
"$RS" dev/rel078-news-merge.R || exit 1
left=$(for f in $(cat "$W/files") $(cat "$LOG"/copied-*.lst); do
  [ -f "$f" ] && grep -l '^<<<<<<< \|^>>>>>>> \|^=======$' "$f"; done)
[ -n "$left" ] && { echo "CONFLICT MARKERS LEFT:"; echo "$left"; exit 1; }

# the working copy's convention is CRLF (core.autocrlf=true, no
# .gitattributes): every text file the merge, the resolver or the
# generators wrote goes back to CRLF, so git sees no line-ending change
# (meme's frmtmb.sample R/draws-brms.R and test-prior-draws.R are LF in
# its tree)
ncr=0
for f in $( { git -C "$G" status --porcelain --untracked-files=all |
              sed -E 's/^.. //'; } | sort -u); do
  [ -f "$f" ] || continue
  case "$f" in
    *.R|*.r|*.md|*.Rmd|*.Rd|*.sh|*.tsv|*.yaml|*.yml|*.jobs|*.ps1|*.cmd|*.patch|*.sed|*.json|*/DESCRIPTION|DESCRIPTION|*/NAMESPACE|NAMESPACE|.Rbuildignore) ;;
    *) continue;;
  esac
  case "$f" in dev/rel078-log/*|dev/ciharden-log/*) continue;; esac
  if ! grep -q $'\r$' "$f"; then
    sed -i 's/$/\r/' "$f"; ncr=$((ncr + 1))
  fi
done
echo "stored with CRLF: $ncr files"
echo "MERGE, RESOLVE AND NEWS DONE"
