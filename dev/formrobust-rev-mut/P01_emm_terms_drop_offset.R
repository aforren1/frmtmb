# The lane's round-0 behavior: emm_terms() strips the offset attribute,
# so emmeans() adds no .offset. column (the blocker).
source("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut/helper.R", local = TRUE)
mut("emm_terms",
    "return(stats::delete.response(tg$targets[[1L]]$lp[[\"terms\"]]))",
    "{ tt <- stats::delete.response(tg$targets[[1L]]$lp[[\"terms\"]]); attr(tt, \"offset\") <- NULL; return(tt) }")
