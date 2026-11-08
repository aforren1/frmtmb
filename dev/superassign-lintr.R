# goodpractice's lint for `<<-`, run alone on core and every extension:
# lintr::undesirable_operator_linter() restricted to the two global
# assignment spellings.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
lnt <- lintr::undesirable_operator_linter(
  op = c("<<-" = "global assignment", "->>" = "global assignment"))
for (d in c("R", file.path(list.dirs("extensions", recursive = FALSE),
                           "R"))) {
  l <- lintr::lint_dir(d, linters = lnt, parse_settings = FALSE)
  cat(sprintf("%s: %d lints\n", d, length(l)))
}
