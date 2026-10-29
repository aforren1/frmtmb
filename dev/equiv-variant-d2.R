# D2: brms's class = "sds" prior against frmtmb's spelling of the same
# prior, class = "sd" with group = the smooth's label (pair p043's data)
spec$prior <- NULL
spec$prior_brms <- "set_prior('student_t(3, 0, 1)', class = 'sds')"
spec$prior_frm <- paste0("set_prior('student_t(3, 0, 1)', class = 'sd', ",
                         "group = 's(x)')")
