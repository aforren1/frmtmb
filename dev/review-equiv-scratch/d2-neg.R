spec$prior <- NULL
spec$prior_brms <- "set_prior('student_t(3, 0, 1)', class = 'sds')"
spec$prior_frm <- "set_prior('student_t(3, 0, 2)', class = 'sd', group = 's(x)')"
