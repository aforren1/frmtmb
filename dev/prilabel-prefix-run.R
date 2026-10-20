# Lane prilabel, fix round 1: run one test file with ONE fix of the
# round taken back out, to see the new test fail before its fix. Each
# mutant puts back the function the lane had before the fix (the text
# of the round-0 lane build, `dev/reviews/2026-10-14-prilabel.md`
# names each defect), in the installed lane namespace.
#   Rscript dev/prilabel-prefix-run.R <fix> <package> <test file>
local({
  a <- commandArgs(TRUE)
  fix <- a[1]
  p <- a[2]
  f <- a[3]
  .libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
              "C:/Users/adf44/source/r/rellib-r9",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
  suppressMessages(library(testthat))
  suppressMessages(library(p, character.only = TRUE))
  ns <- asNamespace("frmtmb")
  put <- function(nm, fn, where = ns) {
    environment(fn) <- where
    utils::assignInNamespace(nm, fn, ns = environmentName(where))
  }
  if (fix == "1") {
    # round 0: one set per factor, refused when the levels differ
    put("student_share_sets", function(re_blocks) {
      st <- which(vapply(re_blocks, is_student_block, NA))
      if (length(st) < 2L) return(re_blocks)
      key <- vapply(re_blocks[st], function(bk) bk[["group_name"]] %||% "",
                    "")
      for (g in unique(key[duplicated(key)])) {
        ids <- st[key == g]
        lv <- re_blocks[[ids[1L]]][["levels"]]
        same <- vapply(re_blocks[ids], function(bk) {
          identical(bk[["levels"]], lv)
        }, NA)
        if (!all(same)) {
          frm_stop("The Student-t terms of grouping factor '", g,
                   "' share one latent scale per level, as in brms, but ",
                   "do not see the same levels", call. = FALSE)
        }
        for (i in ids) {
          re_blocks[[i]][["t_share"]] <- ids
          re_blocks[[i]][["t_level"]] <- seq_along(lv)
          re_blocks[[i]][["t_n"]] <- length(lv)
        }
      }
      re_blocks
    })
  } else if (fix == "2") {
    # round 0: the column name alone, whatever the prefix
    put("sd_coef_shared", function(spec, frame, bk, s) {
      if (!nzchar(s$coef %||% "") || !length(block_sd_idx(bk)) ||
            bk[["covstruct"]] %in% sd_coef_unaddressable) {
        return(character(0))
      }
      if (nzchar(s$group) && !identical(bk[["group_name"]], s$group)) {
        return(character(0))
      }
      if (!all(is.na(block_sd_coefs(frame, bk)))) return(character(0))
      cf <- brms_rename(bk[["cnms"]] %||% character(0))
      cf[cf == s$coef | cf == par_name_bare(s$coef)]
    })
  } else if (fix == "3") {
    # round 0: a smooth's sds_ written inside the block loop
    src <- paste(deparse(get("hyp_env_vals", ns)), collapse = "\n")
    src <- sub("sds_at <- c(sds_at, bi)",
               "put(sds_nm[bi], exp(th[bk[[\"theta_idx\"]][1L]]))",
               src, fixed = TRUE)
    put("hyp_env_vals", eval(parse(text = src)))
  } else if (fix == "4") {
    # round 0: the written order of the terms
    put("re_brms_order", function(re) re)
  } else if (fix == "5") {
    # round 0: every ':' of reformulas's label read as a separator
    put("brms_level_map", function(expr, levels, data, env) NULL)
    put("group_values_brms", function(expr, data, env) {
      gsub(":", "_", as.character(group_values(expr, data, env)),
           fixed = TRUE)
    })
  } else if (fix == "6") {
    # round 0: an interaction's parts in any order
    put("re_group_key", function(g) {
      g <- re_expand_group(g)[[1L]]
      parts <- strsplit(deparse1(g), ":", fixed = TRUE)[[1L]]
      paste(sort(trimws(parts)), collapse = ":")
    })
  } else if (fix == "7") {
    # round 0: no note for a specification that reaches no parameter
    put("unreached_prior_notes", function(fit, user_pl, ri) character(0),
        asNamespace("frmtmb.sample"))
  } else if (fix == "8") {
    # round 0: the sd rows in default_prior()'s order
    put("prior_rows_sd_order", function(tab) tab)
  }
  cat("lib:", p, dirname(find.package(p)), "; mutant of fix", fix, "\n")
  res <- as.data.frame(test_file(f, package = p, env = test_env(p),
                                 reporter = "silent"))
  cat(sprintf("RESULT %s fix=%s pass=%d fail=%d err=%d skip=%d\n",
              basename(f), fix, sum(res$passed), sum(res$failed),
              sum(res$error), sum(res$skipped)))
  bad <- res[res$failed > 0 | res$error, "test"]
  if (length(bad)) cat(paste0("FAILED [", bad, "]"), sep = "\n")
})
