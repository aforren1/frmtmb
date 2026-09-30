source("dev/defects-pre.R")
inh <- brms::inhaler
cat("categorical, cs() as a population term\n")
show(frm(rating ~ treat + cs(period), data = inh, family = categorical()))
cat("categorical, cs() in a group term\n")
show(frm(rating ~ treat + (cs(period) | subject), data = inh,
         family = categorical()))
cat("sratio, cs() in a group term\n")
show(frm(rating ~ treat + (cs(period) | subject), data = inh,
         family = sratio()))
cat("gaussian, cs() in a group term\n")
show(frm(rating ~ treat + (cs(period) | subject), data = inh))
cat("gaussian, cs() population\n")
show(frm(rating ~ treat + cs(period), data = inh))
