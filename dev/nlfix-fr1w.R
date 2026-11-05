source("dev/nlfix-rep.R")
rep1("dev/nlfix-report.R",
r"---(blocks[["a4c"]] <- fence(utils::tail(readLines("dev/nlfix-rev-a4c.log", warn = FALSE), 12)))---",
r"---(blocks[["a4c"]] <- fence(gen_lines("dev/nlfix-rev-a4c.log", "slope|max row")))---")
rep1("dev/nlfix-findings.md",
r"---(  fit and failed the test, so the rule reads the flat direction too. It keeps the
  first fit unless a restart is better.)---",
r"---(  fit and failed the test, so the rule reads the flat direction too. It keeps the
  first fit unless a restart is better, and a restart counts as better only if its
  value is smooth around it (`nl_restart_credible()`): the sweep's beta_binomial row
  read +6.6e12 on a restart that took phi to exp(60) before that check existed.)---")
