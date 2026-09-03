# Lint policy: what the sweep fixes and what it refuses

This record covers the mechanical lint sweep of `R/` done before the
rOpenSci submission. It lists the classes that were swept, the counts
per class per file, and the classes that stay unfixed on purpose.

The scan tool is `lintr` 3.4.0 on R 4.6.1, run per class with
`lint_dir("R", linters = ..., parse_settings = FALSE)`. The `on.exit`
class has no `lintr` linter, so a source scan finds it.

## Scope

Six files were reserved for other work and are not touched by this
sweep: `R/priors.R`, `R/interop.R`, `R/covstruct.R`, `R/objective.R`,
`R/predict.R`, `R/conditional-effects.R`. Their remaining counts are in
the table at the end.

## Classes swept

### 1. `on.exit()` without `add = TRUE`

Each `on.exit()` call without `add = TRUE` replaces every handler that
the function registered before it. The package has six `on.exit()`
calls. Every one of them is the only handler in its function, so no
handler was being replaced and no live bug was found. `add = TRUE` was
added anyway, so that a later handler in the same function cannot
silently delete the first.

| File | Fixed | Sites |
| --- | --- | --- |
| `R/confint.R` | 2 | `hypothesis()`, `plot.frmtmb_hypothesis()` |
| `R/influence.R` | 1 | `plot.frmtmb_influence()` |

The `hypothesis()` site is the one to watch. It arms a shadow-name
note, then calls `UseMethod()`. The handler must survive the dispatch,
and `add = TRUE` does not change that.

### 2. `1:length(x)` and its relatives

Zero hits. `lintr::seq_linter()` reports nothing in `R/`. The only two
`1:` sequences in the package are inside comments, in `R/autocor.R` and
in `R/conditional-effects.R`, where they describe an index range in
prose. The package already uses `seq_len()` and `seq_along()`.

### 3. `sapply()`

Zero hits in package code. The one `sapply(` string in `R/` is in a
roxygen `@examples` block in `R/lca.R`, where the example shows a user
how to summarize a list of refits. Example code is user-facing prose
about a normal R idiom, so it stays.

### 4. `any(is.na(x))` and `any(duplicated(x))`

| File | Class | Fixed |
| --- | --- | --- |
| `R/sandwich.R` | `any(duplicated(...))` | 1 |

`R/sandwich.R:238` becomes `anyDuplicated(pairs[, 1L]) > 0L`.
`anyDuplicated()` returns the integer position of the first duplicate,
or `0L`, so the comparison is needed to keep the `if` condition
logical.

`any(is.na(x))` has zero true hits. `R/ode.R:305` reads
`any(is.na(st) & method != "reset")`, which tests a conjunction and not
a bare `is.na()`, so `anyNA()` does not apply to it.

### 5. Fixed-string patterns

| File | Fixed |
| --- | --- |
| `R/fit.R` | 1 |

`R/fit.R:1440` becomes `grepl("NA/NaN", conditionMessage(w), fixed = TRUE)`.
The pattern is a literal. The slash is not a regular-expression
metacharacter, so the match is identical.

An AST scan of `grep`, `grepl`, `sub`, `gsub`, `regexpr`, `gregexpr`
and `strsplit` finds 39 calls in `R/` that take a literal pattern. 27
of them use metacharacters and keep regular-expression semantics. 11
declare `fixed = TRUE` after this sweep. One candidate is left, in a
reserved file (see the table below).

### 6. `=` used for assignment

Zero hits. `lintr::assignment_linter(operator = c("<-", "<<-"))`
reports nothing. Every assignment in `R/` uses an arrow.

### 7. Long lines

| File | Fixed |
| --- | --- |
| `R/confint.R` | 1 |

`R/confint.R:1594` was 81 columns. The guard clause moves into braces.

## Classes not chased, and why

### Function length and cyclomatic complexity

`goodpractice` reports the objective builder and the frame assembler as
too long and too branched. They stay as they are.

Both functions resolve the model structure in one pass. The frame
assembler reads the formula, the family, the covariance structures and
the data together, because a decision about one of them constrains the
others. The objective builder then writes the tape while those
decisions are still in scope. Splitting either function would move
tape-time invariants across a call boundary, where nothing checks them.
The result would be shorter functions and a weaker guarantee. Length
here is a symptom of a single resolution pass, not of missing
structure.

### `<<-` super-assignment

`lintr::assignment_linter()` reports eight `<<-` operators, in
`R/compat.R`, `R/confint.R`, `R/frame.R`, `R/interop.R` and
`R/priors.R`. Each one is a closure that accumulates rows or names into
a variable in its enclosing function. The pattern is deliberate and
local: the accumulator and the closure are always in the same function
body, and the closure is never returned. Replacing `<<-` with a
`Reduce()` or an explicit list return would change the shape of five
builders for no correctness gain.

### Provably unreachable `1:0`

The package has no `1:n` loop bounds at all, so it has no reachable and
no unreachable zero case. This entry stays in the record so that a
later reviewer does not have to repeat the search.

### Condition-message construction

`goodpractice` prefers conditions that carry a call, and it flags
`call. = FALSE`. `frmtmb` sets `call. = FALSE` on every `stop()` and
`warning()` on purpose.

A frmtmb condition is raised deep inside a builder whose frame the user
never wrote and cannot read. Printing that call adds a line of internal
detail and hides the message. The messages instead name the user-facing
thing that is wrong, such as the formula term, the family, or the
argument. The messages are also load-bearing: `dev/sandwich/count-messages.R`
asserts that all 640 literal `stop()` texts are distinct, and the test
suite matches many of them by text. No sweep may reword them.

### Long lines in `R/compat.R`

208 of the 210 remaining long lines are in `R/compat.R`. They are the
rows of the feature-compatibility registry, and the long part of each
row is the explanation string that the registry returns to the user.
Wrapping those strings would need `paste0()` on every row, which would
make the table harder to read and would risk a whitespace change in
text that users see. The registry keeps one row per line.

### One long line in `R/autocor.R`

`R/autocor.R:708` is 88 columns. It is a message fragment inside a
`stop()` call. Splitting a message literal is a change to load-bearing
text, and the readability gain is nil. It stays.

## Remaining counts in the reserved files

Counts after this sweep. A later round can clear them.

| File | `on.exit` no `add` | seq | `sapply` | `anyNA` | `anyDuplicated` | fixed regex | `=` assign | `<<-` | line > 80 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `R/conditional-effects.R` | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `R/predict.R` | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `R/priors.R` | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 0 |
| `R/interop.R` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 |
| `R/objective.R` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| `R/covstruct.R` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

The three reserved `on.exit()` sites are each the only handler in their
function, so none of them is a live bug. The `R/priors.R` fixed-regex
site is `strsplit(m[3], ",")` at line 103, where the comma is a
literal.

## Verification of this sweep

The sweep touches four files, so the gate is the whole test suite and
not a selection. Every test file ran in its own R session, with
`NOT_CRAN=true` and the worktree source loaded by `pkgload::load_all()`.
93 test files: 5323 pass, 0 fail, 0 error, 6 skip, 0 session failure.

`dev/sandwich/count-messages.R` reports 640 `stop()` literals and 640
distinct texts, the same as before the sweep. A stronger check compares
the message vector of the pristine tree against the message vector of
the working tree and requires them to be identical. That check passes,
so there is no text drift at all.

`roxygen2::roxygenise()` after the sweep leaves `man/`, `NAMESPACE` and
`DESCRIPTION` unchanged, as expected: no roxygen block was edited.

## How to re-run the scan

```r
library(lintr)
setwd("<package root>")
lint_dir("R", linters = list(seq_linter()), parse_settings = FALSE)
lint_dir("R", linters = list(any_duplicated_linter()), parse_settings = FALSE)
lint_dir("R", linters = list(fixed_regex_linter()), parse_settings = FALSE)
lint_dir("R", linters = list(assignment_linter(operator = c("<-", "<<-"))),
         parse_settings = FALSE)
lint_dir("R", linters = list(line_length_linter(80L)), parse_settings = FALSE)
```

`lintr` has no `on.exit` linter. Find that class with a source scan for
`on.exit(` lines that do not also contain `add =`.
