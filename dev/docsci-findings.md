# Lane wt-docsci: the documentation site moves to CI

Lane: `wt-docsci`. Worktree: `C:\Users\adf44\source\r\frmtmb-wt-docsci`.
Dates: 2026-09-28 into 09-29. Nothing here is committed.

## What the problem was

The site at <https://aforren1.github.io/frmtmb/> was built by
`dev/release/run-docs.ps1` on one Windows box and committed to `docs/`.
Three things follow from that, and all three were measured rather than
assumed:

- **The published site is a function of one machine's library.**
  `run-docs.ps1` hard-codes `R_LIBS` to `C:/Users/adf44/source/r/rellib-r3`
  and a Windows pandoc path. Nobody else can reproduce the site.
- **Every rebuild is a large commit.** `git ls-files docs | wc -l` gives
  1,304 files and `du -sh docs` gives 41 MB.
- **Nothing deletes what a rebuild stops producing,** and one thing in
  the served tree is not produced by any build at all. See **What a
  clean rebuild loses**. That section is the real work of this lane and
  it carries the one item that blocks the cutover.

Pages was and still is on the legacy builder. Read on 2026-09-28:

```
$ gh api repos/aforren1/frmtmb/pages
{"status":"built","build_type":"legacy",
 "source":{"branch":"main","path":"/docs"}, ...}
```

## What I built

`dev/release/build-docs.R`
: the build, the layout and every guard. One R file, no machine
  assumptions.

`.github/workflows/pkgdown.yaml`
: builds the eight sites on `ubuntu-latest` and deploys to Pages. The
  deploy job is conditional on the cutover.

`dev/release/run-docs.ps1`
: rewritten as a thin wrapper: it sets this box's paths and calls
  `build-docs.R`. Kept rather than deleted, because the Windows library
  and pandoc paths in its comments are real knowledge.

`CONTRIBUTING.md`
: a new **Documentation** section: how to build locally, what the guards
  check, which articles do not run their code and why.

`.gitignore`
: the scratch build destinations only. A `docs/` line was withdrawn at
  consolidation, 2026-09-29: while `docs/` is still tracked, the rule does
  not unignore the files a rebuild ADDS, so a docs commit would silently
  leave every new reference page out. It is added in step 6 instead.

`extensions/*/_pkgdown.yml`
: seven comments said "Build it with `dev/build-docs.R`", a file that did
  not exist. They now name `dev/release/build-docs.R`, and the paragraph
  is rewrapped to 80 columns around the longer path. Comments only, and
  the `destination:` of all seven is byte for byte what it was.

`.Rbuildignore` needs no change: it already carries `^docs$` and `^dev$`.

Scratch, all prefixed `docsci-` and all under `dev/`. They are the
record of how each number below was produced, which is why they are
left in the tree rather than deleted.

| file | what it did |
|---|---|
| `docsci-verify.R` | generated the per-site block in this document |
| `docsci-validate-all.ps1` | the full local run: build, then guards |
| `docsci-run-local.ps1` | the eight-site build on this box |
| `docsci-guard-tests.ps1` | the eight guard cases |
| `docsci-dryrun-deps.R` | the extension dependency parser, no installs |
| `docsci-dryrun-core-suggests.R` | the core Suggests parser, no installs |
| `docsci-dryrun-gates.R` | the derived article gate list |
| `docsci-check-articles.R` | `articles:` blocks against the `.Rmd` files |
| `docsci-run-ci-test.R` | `test-ci-siblings.R`, one file one process |
| `docsci-fix-yml-comment.R` | rewrapped the one comment paragraph, safely |
| `docsci-inspect-redirect2.R`, `-sitemap.R` | read the pkgdown internals |
| `docsci-compare-search.R` | the two `search.json` files, entry by entry |
| `docsci-search-delta.R` | the final index and sitemap delta, generated |
| `docsci-verify-final.txt` | the first generated block, as generated |
| `docsci-search-delta.txt` | the second generated block, as generated |
| `docsci-check-redirects.R` | classified the 105 pages a clean build drops |
| `docsci-gen-redirects.R` | wrote `docsci-core-redirects.yml` |
| `docsci-core-redirects.yml` | the 104 redirects, for the owner to decide |
| `docsci-build.log` | the log of the run this document cites |
| `docsci-build-run1.log` to `-run3.log` | earlier runs, for comparison |
| `docsci-validate-all.out` | the guard transcript, every line AS EXPECTED |

### The layout is unchanged

Each extension's `_pkgdown.yml` carries `destination: ../../docs/<pkg>`.
`build-docs.R` reads that file and takes the LAST SEGMENT of
`destination` as the subdirectory, so the tree it writes is the tree that
is published today, and the source of truth stays in the same file it
was in. It refuses a `destination` whose last segment is not the package
name, because that would put a site where the navbar does not link.

### Stan and brms do not run in CI, and the site is not degraded

This is the question the task said not to answer by assumption, so it
was answered against the committed `docs/` tree.

| article | gate | `^#>` lines in `docs/` |
|---|---|---|
| `brms-migration.Rmd` | whole-document `eval = FALSE` | 0 |
| `sample/brms-posterior.Rmd` | whole-document `eval = FALSE` | 0 |
| `sample/posterior-diagnostics.Rmd` | whole-document `eval = FALSE` | 0 |
| `sample/sampling.Rmd` | per-chunk `eval = FALSE` | 20, sampler-free only |

`grep -c "^#>" docs/articles/brms-migration.md` is 0. The two
`FRMTMB_BRMS_FIT_TESTS` mentions, in
`vignettes/bayesian-cognitive-modeling.Rmd` and
`vignettes/reinforcement-learning.Rmd`, are inside fenced prose blocks
that describe the TEST suite. Neither is a chunk, and neither is
evaluated. No `_pkgdown.yml` names an article that is not a vignette,
and there is no `pkgdown/` directory in any of the eight packages.

So: nothing on the published site depends on Stan, brms or tmbstan, and
the workflow installs no Stan toolchain and sets no `FRMTMB_*` gate. The
workflow header says this, and `CONTRIBUTING.md` says what to do if an
article ever changes: reuse the Stan setup in `brms-likelihood.yaml`
rather than let the article render short.

Every article on the site is a vignette. `dev/docsci-check-articles.R`
compared each `articles:` block against the `.Rmd` files: the core names
10 and has 10, and the seven extensions name none, so pkgdown takes all
11 of theirs by default. 21 articles, 21 files, no name without a file
and no file without a page. There is no `pkgdown/` directory in any of
the eight packages, so there is no article that is not a vignette and no
`eval` gate hiding anywhere else.

`sampling.Rmd` already carries a visible note about why its chunks do
not run, so nothing had to be added. Quoting it:

> Most chunks on this page are shown and not run, and speed is not the
> reason ... The reason is that tmbstan and rstan are `Suggests`, so a
> build machine need not have them, and that a tmbstan built against
> the wrong StanHeaders samples the wrong density in silence.

### What WOULD degrade the site, and the guard for it

Every article gates its fitted output on `requireNamespace()`. A build
machine missing one of those packages renders the page with its numbers
dropped and exits 0. That is the silent regression, and it is the one
the old driver had no defense against.

`build-docs.R` reads the gate names out of the article sources with a
regex over `requireNamespace("...")`, so the list cannot go stale, and
`--require-articles` turns a missing one into a failure BEFORE the build
starts. The workflow always passes it. A local build gets a warning
instead, so a quick local rebuild still works.

On this machine the script found 17 gate packages and none missing.
`RTMBode` is the only one that is not on CRAN; the workflow adds
`extra-repositories: 'https://kaskr.r-universe.dev'`, copying
`check-frmtmb-ode.yaml`. Without it the whole `frmtmb.ode` article goes
quiet.

## The guards, and each one seen to fire

Before anything is built, in this order, so that a refusal leaves the
site already in the destination intact:

- **P1** every package installed, and its installed version equal to its
  `DESCRIPTION`. pkgdown renders the INSTALLED package, so a stale
  install publishes stale output and says nothing.
- **P2** every `destination` in an extension `_pkgdown.yml` ends in the
  package name.
- **P3** with `--require-articles`, every derived article gate is
  installed.
- **P4** the destination is either empty or already a pkgdown site
  before it is emptied.

After the build, on files the run produced and not on the loop that
produced them:

- **G1** every site exited 0.
- **G2** every site has a non-empty `index.html` AND a
  `reference/index.html`. The second catches a site that died after its
  home page, which is exactly the shape `docs/frmtmb.ddm/` has.
- **G3** the number of sites passing G2 is at least `--min-sites`,
  default 8.
- **G4** the root reindex step exited 0, and `search.json` exists at the
  root. Full builds only.
- **G5** the root `search.json` holds at least one path under EVERY
  extension subdirectory, and the root `sitemap.xml` holds at least one
  URL under every one. Full builds only. This is the guard for the
  regression in **Difference 1** below, and it exists because the
  existence check in G4 passed happily while the root search index had
  lost 978 of its 2,088 entries.
- **G6** every `href` in the core navbar Extensions menu resolves to a
  built subdirectory. Full builds only.
- **G7** stray directories in the destination are reported.

G4 to G6 are about the whole tree, so a deliberately partial build would
fail them for the wrong reason and `--pkgs` would be useless for a quick
local rebuild. They run on a full build and on demand with
`--check-links`. A partial build prints, instead, that the whole-tree
checks were skipped and that the tree must not be deployed.

The workflow adds a second count in bash, off the filesystem, with the
expected number derived from `ls extensions/*/_pkgdown.yml` and an
absolute floor of 8. Neither side of that comparison is the R script's
own tally. As a check on the check: run that same `find` against the
committed `docs/` tree and it returns 9, not 8, because of
`docs/frmtmb.ddm/`. The count can tell the two trees apart.

## What a clean rebuild loses, which is the real work of this lane

The first eight-site build looked fine: 8 of 8, every navbar link live,
no article with fewer output lines than the published one. Then the
trees were compared file by file, and two differences turned up that no
count would have shown. A retraction is recorded below as well, because
the first reading of one of them was wrong.

### Difference 1: the unified search box and sitemap at the root

**FIXED, on the second attempt. The first fix failed loudly, which is
recorded here rather than quietly replaced.**

pkgdown's `build_search()` indexes the HTML it finds in the DESTINATION,
not only the package it is building, and `pkgdown:::build_sitemap()`
walks the destination the same way. The old driver built core first, into
a `docs/` that still held the PREVIOUS run's subsites, so the root
`search.json` and `sitemap.xml` picked all of them up. Measured on the
two search indexes with `jsonlite`:

| | committed `docs/` | core-first clean build |
|---|---|---|
| entries | 2,088 | 1,110 |
| paths only on this side | 125 | 1 |
| characters | 2,015,836 | 1,101,555 |

All 125 are subsite pages: 29 under `frmtmb.eam/reference`, 20 under
`frmtmb.sample/reference`, 14 under `frmtmb.learn/reference`, and so on
down to the subsite article pages, worth 852,999 of the committed index's
characters. The single path only in the new index is
`reference/log_lik.html`, a core page the committed index misses. The
root `sitemap.xml` says the same thing: 497 `<loc>` entries, of which 268
are subsite URLs.

So the root search box on the live site finds extension pages, by
accident, and a clean build with core first silently takes that away.

**The obvious fix does not work.** Building core LAST was tried first. It
failed at once, on all eight sites' worth of guards:

```
===== SITE frmtmb -> .../dev/docsci-site =====
Error in `check_dest_is_pkgdown()`:
! '.../dev/docsci-site' is non-empty and not built by pkgdown
! ... use `pkgdown::clean_site(force = TRUE)` to delete its contents.
===== FAILED exit 1 2s =====
DOCS built 7 of 8 (minimum 8)
DOCS FAILED
  - site frmtmb exited 1
  - no usable index.html for frmtmb at .../dev/docsci-site/index.html
  - built 7 sites, expected at least 8
  - no search.json at the root of .../dev/docsci-site
```

`pkgdown:::check_dest_is_pkgdown()` refuses a non-empty destination whose
ROOT has no `pkgdown.yml`. Seven subdirectories are not a pkgdown site.
Note what happened: four separate guards fired and named the cause. That
is the guard design working, and it is why the order change did not ship
as an unverified one-liner.

**The fix that does work** is core first, the extensions after it, and
then the root index regenerated over the finished tree:

```r
p <- pkgdown::as_pkgdown(core, override = list(destination = dest))
pkgdown::build_search(p)          # exported
pkgdown:::build_sitemap(p)        # not exported
```

`build_search()` is exported. `build_sitemap()` is not, so the script
looks it up with `get0()` and fails with a message naming the problem if
a future pkgdown renames it, rather than skipping the sitemap in silence.
The step runs as its own child process, like every site build, and only
on a full build: reindexing a partial tree would write a root index for a
site that is not all there.

**And a guard, because the loss was silent.** G5 reads the finished
`search.json` and `sitemap.xml` and requires a path under every extension
subdirectory in each. The existence check that was already there, G4,
passed the whole time the index was missing 978 of its 2,088 entries. A
guard on the file being present is not a guard on the file being right.

### Difference 2: 104 redirect pages that nothing regenerates

**NOT FIXED. This blocks the cutover.**

`docs/reference/` holds 413 HTML pages and a clean core build produces
308. The difference is one way: 105 pages only in the committed tree, 0
only in the new one, and every extension's own count matches exactly
(coupling 9, eam 38, latent 10, learn 15, ode 7, sample 134, spline 8).

Of those 105, **104 are meta-refresh redirect pages** and one,
`frm_family.html`, is a real page for a topic core no longer has.

The 104 are not accidents. `git log` on one of them lands on
9b9011a4, 2026-09-05, "retire orphaned core reference pages": the owner
replaced 518-line pages with 7-line redirects when the functions moved
into the extensions. They are hand-written files. pkgdown does not
produce them and will not reproduce them, so the clean at the start of
the new build deletes all 104 and those URLs begin answering 404.

pkgdown's own `redirects:` key writes the same page from config. Reading
`pkgdown:::build_redirect` rather than recalling it, the target is
formed as `sprintf("%s/%s%s", url, prefix, new)` with `prefix` empty for
a released site, so `new` is relative to the core site root and a
subsite target spells itself `frmtmb.sample/reference/x.html`. Every one
of the 104 can be expressed that way. `dev/docsci-gen-redirects.R`
generates the block into `dev/docsci-core-redirects.yml`:

- **53** have a target that a clean build produces. Those lines can be
  pasted into the core `_pkgdown.yml` as they stand.
- **51** point at a page that is itself one of the 104, for example
  `reference/draws-diagnostics.html` and
  `reference/frmtmb-loo-refusals.html`. Those chains already end
  nowhere once the tree is cleaned; their real target is a
  `frmtmb.sample` page, and choosing it is a decision per name rather
  than a mechanical rewrite. They are in the file commented out.

I did not paste either group into `_pkgdown.yml`. 51 of them need a
judgement about where an old URL should land, and inventing 51 answers
silently is worse than handing over the list. This is step 2 of the
runbook for that reason.

### Retracted: the claim that those pages are in the search index

My first reading was that the 104 redirect pages were being served in
the core search index. It came from
`grep -o "reference/frm_sample.html" docs/search.json`, which returns 2,
and the reading was wrong: the paths in `search.json` are absolute URLs,
so that pattern also matches
`.../frmtmb.sample/reference/frm_sample.html`. Anchored on the core
prefix instead, exactly **one** of the 105 is in the core index, and it
is `frm_family.html`, the one that is a real page rather than a
redirect. The redirect stubs carry `meta name="robots" content="noindex"`
and almost no text, and they are not indexed. The 1 MB of index that the
new build drops is Difference 1 and nothing to do with these.

### Two more kinds of stale, both harmless and both removed by the clean

- `docs/frmtmb.ddm/`, the site of the package renamed to `frmtmb.eam`,
  last touched 2026-09-05. It is a husk: `index.html` exists,
  `reference/index.html` does not, no articles. It answers 200 today.
- `docs/articles/ode_files/`, in the CORE site, holding the same three
  PNG files as `docs/frmtmb.ode/articles/ode_files/`. The core package
  has no `ode.Rmd`.

## Every guard, seen to fire, with a control that must not

`dev/docsci-guard-tests.ps1` constructs eight cases. Each case that must
FAIL is paired with one that must PASS, because a check that fires on
correct input is worse than no check. The harness prints the expected
outcome beside the observed one and labels each `AS EXPECTED` or `WRONG`,
so the transcript in `dev/docsci-validate-all.out` cannot be read
optimistically.

**A**, all 8 sites, `--min-sites 8`: must pass. Nothing fired, `DOCS ok`.

**C**, 1 site, `--min-sites 1`: must pass. Nothing fired. The control for
B and D, so that their failure is the count and not the partialness.

**D**, 1 site, `--min-sites 8`: must fail. G3 alone, and it said
`built 1 sites, expected at least 8`.

**B**, the 7 extensions, `--min-sites 8`, the case the task named: must
fail. G3 alone, `built 7 sites, expected at least 8`.

**E**, a destination holding `precious.txt`: must fail. P4 fired and the
file survived.

**F**, a planted unsatisfiable gate with `--require-articles`: must fail.
P3 fired, and 0 entries were left in the destination, which is the proof
that it refused before the clean and before any build.

**G**, the same gate without `--require-articles`: must pass. A warning,
then a clean build. The control for F.

**H**, `destination:` bent to `frmtmb.coupled`: must fail. P2 fired
during discovery, before any site was built.

**I**, a navbar `href` bent to `frmtmb.coupled`: must fail. G6 fired
alone, naming the href and the directory it wants:
`navbar href https://aforren1.github.io/frmtmb/frmtmb.coupled/ has no
built site at frmtmb.coupled/index.html`.

The transcript is `dev/docsci-validate-all.out` and every line of it
ends in `AS EXPECTED`. H and I bend a `_pkgdown.yml` and restore it, and
F and G plant a vignette and remove it; the harness asserts each
restoration and all four reported true. `git diff` on the core
`_pkgdown.yml` is empty, and on each extension's it is the comment
paragraph this lane rewrote and nothing else.

**Two of the guards were wrong on their first try, and both are recorded
here rather than quietly fixed.** Neither was found by review.

1. **P4 failed CLOSED.** The clean refusal treated "no `pkgdown.yml` or
   `index.html` at the root" as "not a pkgdown site", which is false for
   a destination holding only extension subdirectories. That is exactly
   what a partial build of the extensions produces, so case D refused to
   clean and failed on the refusal instead of on the count arm it was
   written to test. The case that was supposed to test the count guard
   was passing for the wrong reason. Fixed by also accepting a top level
   whose every entry is a directory containing its own `pkgdown.yml`;
   case E proves the refusal still works on a plain file.

2. **Case G was a bad TEST, not a bad guard.** The probe vignette was
   planted in the extension being built, so pkgdown tried to render a
   file with no `VignetteIndexEntry` and the site build died. G reported
   `WRONG`. Moving the probe into the CORE vignettes directory, while the
   case builds only `frmtmb.coupling`, means the probe is SCANNED by the
   gate detector and never rendered, which is what the case was for.

A third failure, the core-last build, is in **Difference 1** above. Three
guard or test defects in one lane, all caught by running the case rather
than by reading the code, is the standing rule about constructing the
absent case earning its place again.

## The build, measured

Generated by `dev/docsci-verify.R` from the tree the final run produced,
and pasted verbatim. `dev/docsci-build.log` is that run's log.

<!-- generated by dev/docsci-verify.R; do not edit by hand -->
```
site C:/Users/adf44/source/r/frmtmb-wt-docsci/dev/docsci-site
site                  index  refpg  artcl outlines
frmtmb (root)           yes    yes     10     580
frmtmb.coupling         yes    yes      1      46
frmtmb.eam              yes    yes      2     252
frmtmb.latent           yes    yes      1      90
frmtmb.learn            yes    yes      1      35
frmtmb.ode              yes    yes      1      73
frmtmb.sample           yes    yes      3      20
frmtmb.spline           yes    yes      2     130
site roots with index.html and reference/index.html: 8
search.json at the root: yes (2151574 bytes)
reference pages, all sites: 529
files, all sites: 1166

core navbar Extensions menu against the built tree:
  frmtmb.coupling    built
  frmtmb.eam         built
  frmtmb.latent      built
  frmtmb.learn       built
  frmtmb.ode         built
  frmtmb.sample      built
  frmtmb.spline      built

rendered output lines, this build against .../docs:
  article                                            ref     new
  frmtmb.eam/articles/ddm.md                         202     206
articles with FEWER output lines than the reference: 0
```

The run's own tail, from `dev/docsci-build.log`:

```
DOCS seconds per site: frmtmb 354, frmtmb.eam 67, frmtmb.latent 21,
  frmtmb.ode 33, frmtmb.spline 23, frmtmb.coupling 17,
  frmtmb.sample 39, frmtmb.learn 27
DOCS total 626s (10.4 min)
DOCS built 8 of 8 (minimum 8)
DOCS ok
```

Ten and a half minutes, with four other R processes on the machine, and
the core site is 354 s of it. The workflow's 90 minute timeout is mostly
for the dependency install, which this box does not measure.

**Every article renders at least as much as the published one.** Only one
article differs at all, and it went UP: `frmtmb.eam/articles/ddm.md` has
206 output lines against the published 202. Zero articles have fewer.
That is the answer to the question about a CI build silently rendering a
gated-off version of an article.

### The root index, after the reindex fix

Generated by `dev/docsci-search-delta.R`, pasted verbatim:

```
search entries: committed 2088, new 2091
subsite paths: committed 899, new 900
only in committed:
  https://aforren1.github.io/frmtmb/reference/frm_family.html

only in new:
  https://aforren1.github.io/frmtmb/reference/log_lik.html

sitemap locs: committed 497, new 391
subsite locs: committed 268, new 267
only in committed: 106
of those, under the CORE reference: 105
first 5:
  https://aforren1.github.io/frmtmb/frmtmb.ddm/index.html
  https://aforren1.github.io/frmtmb/reference/as.array.frmtmb_draws.html
  https://aforren1.github.io/frmtmb/reference/as.mcmc.frmtmb_draws.html
  https://aforren1.github.io/frmtmb/reference/as.mcmc.html
  https://aforren1.github.io/frmtmb/reference/as_draws.frmtmb_draws.html

only in new: 0
```

The search index is whole: 2,091 entries against 2,088, 900 subsite paths
against 899, and the whole difference is one path each way.
`reference/frm_family.html` is only in the committed index because that
stale page is only in the committed tree. `reference/log_lik.html` is
only in the new index, a core page the published one misses.

The sitemap's 106 missing entries are Difference 2 and the stale
`frmtmb.ddm` site, and nothing else: 105 of them are core reference pages
a clean build does not produce, which is the 104 redirect stubs plus
`frm_family.html`, and the 106th is `frmtmb.ddm/index.html`. Nothing is
in the new sitemap that is not in the old one. Deciding the 104 redirects
closes this gap along with the 404s.

## The cutover runbook, for the repository owner

None of this was done by the lane. Step 2 is a judgement and step 3 is
the owner's setting, and the workflow is written so that it does no harm
before either.

1. **Merge this branch first.** The workflow has to be on `main` before
   `workflow_dispatch` can offer it. Merging changes nothing about what
   is served: Pages is still `legacy`, so the `deploy` job skips itself
   and prints a warning naming this file. The `build` job runs on the
   merge commit and is the first proof that the site builds on a runner.

2. **Decide the 104 redirects, BEFORE anything is deployed.** This is
   the one step that is a judgement and not a command. Read
   `dev/docsci-core-redirects.yml`.

   - Paste the 53 uncommented entries into the core `_pkgdown.yml` under
     a top-level `redirects:` key. They are correct as generated.
   - For the 51 commented ones, decide where each old URL should land.
     Their current target is another page that a clean build does not
     produce, so the honest answer for most of them is the
     `frmtmb.sample` page the function actually lives on now.
   - Or decide that these URLs may 404. That is a defensible choice and
     it is faster; it is just not a choice a lane should make silently.

   Then run the build locally and check the count came back:

   ```
   Rscript dev/release/build-docs.R --dest=/tmp/site
   ls /tmp/site/reference/*.html | wc -l   # 308 before, 412 after
   ```

3. **Switch Pages to the workflow builder.**

   ```
   gh api -X PUT repos/aforren1/frmtmb/pages -f build_type=workflow
   gh api repos/aforren1/frmtmb/pages --jq '.build_type'   # workflow
   ```

   Or Settings, Pages, Build and deployment, Source, GitHub Actions.
   The moment this flips, the legacy `pages build and deployment` job
   stops running and NOTHING serves the site until step 4 finishes. The
   gap is one build, about 30 to 60 minutes on a cold runner. Do it when
   a gap is acceptable.

4. **Run the workflow by hand and watch it.**

   ```
   gh workflow run pkgdown.yaml
   gh run watch "$(gh run list --workflow pkgdown.yaml \
                     --limit 1 --json databaseId --jq '.[0].databaseId')"
   ```

   Read the `Count the built sites again, independently` step. It prints
   `site roots expected 8, found 8` and the list of directories. If it
   prints anything else, stop: the deploy will not have run.

5. **Check the deployed site.** The four things that would break, and
   what to look at:

   - the Extensions dropdown in the navbar, all seven entries, from
     <https://aforren1.github.io/frmtmb/>. The build already asserts
     each href has a directory, so this is checking the deploy, not the
     build.
   - <https://aforren1.github.io/frmtmb/frmtmb.ddm/> returns 404. It is
     200 today. This is the stale directory going away, and it is the
     single clearest signal that the deployed tree is the built tree and
     not the committed one.
   - search, BOTH KINDS. In the root site's search box type
     `frm_bootstrap`, a core topic, and then `frm_sample`, which lives in
     `frmtmb.sample`. Both must return a hit, and the second must land on
     a `frmtmb.sample/reference/` URL. Cross-package search at the root
     is the thing **Difference 1** nearly lost; G5 asserts the index
     holds paths under every extension, but only the deployed site proves
     the browser can read a 2 MB index.
   - one article with fitted output, for example
     <https://aforren1.github.io/frmtmb/frmtmb.ode/articles/ode.html>.
     It must show `#>` output blocks. An empty one means `RTMBode` did
     not install, which `--require-articles` should have failed on
     first.

6. **Only then, remove `docs/` from the repository.**

   ```
   git rm -r --cached docs
   echo "docs/" >> .gitignore
   git add .gitignore
   git commit -m "docs: stop tracking the built site"
   ```

   `.Rbuildignore` already has `^docs$`. The `.gitignore` line is added
   HERE and not earlier: while `docs/` is tracked it would hide every
   file a rebuild adds from `git add`. That the ignore rule bites once
   the file is out of the index was checked rather than assumed, on the
   branch before the line was moved here:

   ```
   $ git check-ignore -v --no-index docs/index.html
   .gitignore:32:docs/    docs/index.html
   ```

   Keep the local directory or delete it; git no longer cares. Doing
   this BEFORE step 5 would leave no way back: the committed tree is
   the rollback.

7. **Rollback, if step 5 fails.**

   ```
   gh api -X PUT repos/aforren1/frmtmb/pages \
     -f build_type=legacy -f 'source[branch]=main' -f 'source[path]=/docs'
   ```

   This works only while `docs/` is still tracked, which is why step 6
   is last.

## Open questions for the user

1. **`docs/frmtmb.ddm/` is served today and has been since 2026-09-05.**
   Does anything link to it? If an external page does, a redirect page
   is a two-line addition to the build script. I did not add one,
   because I have no evidence anyone follows it.

2. **The 104 redirects are step 2 of the runbook and the only thing
   that blocks it.** The question for you is the 51 whose target is
   itself a page a clean build does not produce. Should each go to its
   `frmtmb.sample` page, or should those URLs be allowed to 404? I have
   an opinion, which is that a URL the owner deliberately preserved
   three weeks ago should keep working, but it is your call and the list
   is generated and waiting.

   A related question I did not answer: whether pkgdown could be made to
   do this without config at all, by putting the 104 files under
   `pkgdown/assets/reference/`, which `pkgdown:::copy_assets` copies to
   the site root on every build. That works and needs no decision about
   targets, at the cost of 104 small files in the repository instead of
   104 lines of YAML. `redirects:` is the better shape; the asset route
   is the one that needs no thinking. Say which you want.

3. **Timeout.** The workflow is set to 90 minutes. The local eight-site
   build is the measurement behind that; see the generated block above
   for what it took here. A first CI run also pays for the dependency
   install, which is large: the full `Suggests` of eight packages. If
   the first run comes close to 90 minutes, the next thing to add is an
   `actions/cache` step on the R library, the way
   `brms-likelihood.yaml` caches its compiled models.

4. **Should the workflow also run on pull requests?** It does not today.
   A PR build would catch a broken article before merge, at the cost of
   a 30 to 60 minute job on every PR touching `R/`. I left it off
   because `R-CMD-check` already builds the vignettes on four platforms,
   which is where a broken vignette surfaces first. Say the word and it
   is one `pull_request:` block, with the deploy job already gated.

5. **RESOLVED, not an open question. A `drmTMB` reading that changed
   under me, and why.** Early in the session a `setdiff` against
   `rownames(installed.packages())` with `R_LIBS` set to
   `rellib-r3;<user library>` reported `drmTMB` missing. About forty
   minutes later the same query with the same `R_LIBS` reported it
   present at 0.7.0 in `C:/Users/adf44/source/r/rellib-r3`. I recorded
   that as a possible hollow-directory event and asked for the library to
   be looked at.

   **Do not open a library-loss investigation. The cause is known and
   benign.** The consolidating session installed `drmTMB` 0.7.0 into
   `rellib-r3` at about 22:55 on 2026-09-28, because `R CMD check` on
   core stopped at "Package suggested but not available: 'drmTMB'"
   (`dev/release/check.log` line 33). `rellib-r3` is that session's own
   library, not a lane's, so the write was deliberate and the two
   readings are simply before and after it. Both readings were correct.

   Nothing else changes. It has no bearing on the site either way: no
   article and no `man/` page mentions `drmTMB`, and it is not in the
   derived gate list.

6. **A release deploys the TAGGED tree.** The `release: published`
   trigger is in the workflow because the task asked for it, but a
   release cut from an old tag would roll the published documentation
   back. The push trigger already covers a release made from the tip of
   main. Say the word and the two lines come out.

## What I did not do, and why

- **I did not remove `docs/` and I did not touch the Pages setting.**
  Both are the owner's, and both are one-way without the other.
- **I did not delete `run-docs.ps1`.** It is now 45 lines that set this
  machine's paths and call the R script. Deleting it would throw away
  the `R_MAKEVARS_USER` and pandoc knowledge in its comments for no
  gain.
- **I did not add a Stan toolchain to the workflow.** Measured: no
  article renders any Stan or brms output. Adding one would cost tens
  of minutes per run and change nothing on the page.
- **I did not lint the workflow with `actionlint`.** It is not on PATH
  on this machine, and neither is `npx`, `node`, `yamllint` or a real
  `python3` (`python3` resolves to the Windows Store stub). The YAML was
  checked with `yaml::read_yaml()` and with `yaml::yaml.load()` reading
  the `on:` block key by key, which is reported below. That is weaker
  than `actionlint` and I am saying so rather than implying a clean
  lint.
