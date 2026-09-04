# Forced-response habit data from Hardwick et al. (2019)

Trial-level records from the timed-response assessment of Hardwick,
Forrence, Krakauer and Haith (2019). Participants learned a mapping from
four symbols to four keys, then learned a revised mapping in which two
of the four symbols changed key. In the assessment they had to respond
at a moment fixed by the experimenter, which sets the time available to
prepare a response. The interest is in what happens at short preparation
times, when the original mapping can win over the revised one.

## Usage

``` r
habit_prep
```

## Format

A data frame with 28,940 rows and 10 columns:

- participant:

  Factor. Person identifier, `e1-NN` for experiment 1 and `e2-NN` for
  experiment 2. Experiment 1 identifiers occur in both the `minimal` and
  the `4day` group.

- group:

  Factor with levels `minimal`, `4day`, `20day`. Amount of practice on
  the original mapping before the remapping.

- subject_code:

  Integer. The participant directory in the source deposit, which
  numbers the same person differently in each group.

- trial:

  Integer. Position within the timed-response assessment. Some indices
  are absent because the corresponding trial was not recorded.

- prep_time:

  Numeric. Preparation time in seconds, from stimulus onset to the
  required response. Values are negative when the response came before
  the stimulus appeared.

- remapped:

  Logical. `TRUE` for the two stimuli whose key changed between the
  original and the revised mapping, `FALSE` for the two that kept their
  key.

- stimulus:

  Integer. Symbol shown. The two experiments use different symbol sets,
  so the codes are not comparable between them.

- target_key:

  Integer in 1:4. Correct key under the revised mapping.

- response_key:

  Integer. Key pressed. A few trials record a value outside 1:4, which
  is a hardware read failure rather than a key.

- response:

  Factor with levels `correct`, `habit`, `other`. The response category
  the model is fit to: `correct` is the revised mapping, `habit` is the
  original mapping, and `other` is either of the two remaining keys.
  `NA` where `response_key` lies outside 1:4.

## Source

Hardwick RM, Forrence AD, Krakauer JW, Haith AM (2019). Time-dependent
competition between goal-directed and habitual response preparation.
*Nature Human Behaviour*, 3, 1252-1262.
[doi:10.1038/s41562-019-0725-0](https://doi.org/10.1038/s41562-019-0725-0)

Derived from the authors' data deposit, Hardwick RM and Haith AM (2019),
"Time-Dependent Competition Between Goal-Directed and Habitual Response
Selection", <https://osf.io/3fjez/>. Used under the MIT License,
copyright (c) 2019 Robert Michael Hardwick and Adrian Mark Haith. See
`system.file("COPYRIGHTS", package = "frmtmb")`.

`data-raw/habit_prep.R` in the package sources rebuilds this object from
the deposit.

## Details

Three groups differ in how much the original mapping was practiced
before the remapping: `minimal` is experiment 1 before training, `4day`
is the same people after four days of practice, and `20day` is a
separate cohort trained for twenty days. Because `minimal` and `4day`
are the same people measured twice, `participant` repeats across those
two groups.

Each trial is assigned to the remapped or the unchanged stimulus set by
its trial index.

Only the two remapped stimuli carry a habitual response option, so the
model in
[`vignette("habit")`](https://aforren1.github.io/frmtmb/articles/habit.md)
is fit to `subset(habit_prep, remapped)`.

## Examples

``` r
str(habit_prep)
#> 'data.frame':    28940 obs. of  10 variables:
#>  $ participant : Factor w/ 36 levels "e1-01","e1-02",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ group       : Factor w/ 3 levels "minimal","4day",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ subject_code: int  101 101 101 101 101 101 101 101 101 101 ...
#>  $ trial       : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ prep_time   : num  1.163 1.178 0.762 0.136 0.823 ...
#>  $ remapped    : logi  TRUE TRUE TRUE TRUE FALSE FALSE ...
#>  $ stimulus    : int  3 4 3 3 1 2 1 2 4 2 ...
#>  $ target_key  : int  2 4 2 2 3 1 3 1 4 1 ...
#>  $ response_key: int  2 4 2 1 3 1 4 2 1 1 ...
#>  $ response    : Factor w/ 3 levels "correct","habit",..: 1 1 1 3 1 1 3 3 3 1 ...
with(subset(habit_prep, remapped), table(group, response))
#>          response
#> group     correct habit other
#>   minimal    3777   616  1094
#>   4day       3379  1125   981
#>   20day      2149   776   570
```
