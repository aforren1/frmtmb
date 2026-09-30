for (f in c("extensions/frmtmb.sample/man/frm_sample.Rd","extensions/frmtmb.sample/man/posterior_epred.Rd","man/variables.Rd")) { cat(f, "checkRd problems:", length(tools::checkRd(f)), "\n") }
