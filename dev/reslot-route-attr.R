# Lane reslot, fix round 2: does default_prior()'s route attribute
# survive the edits a user makes to the table? Usage:
# Rscript dev/reslot-route-attr.R before|after. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
dd <- re_data(1L)
f <- y ~ x + (1 + x | g)
p <- default_prior(f, data = dd)
show <- function(label, x) {
  r <- tryCatch(x, error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-36s ERROR %s\n", label, conditionMessage(r)))
  } else {
    cat(sprintf("%-36s class %-40s route %s\n", label,
                paste(class(r), collapse = "/"),
                format(attr(r, "route") %||% "<none>")))
  }
}
show("default_prior()", p)
q <- p; q$prior[q$class == "b" & q$coef == "x"] <- "normal(0, 1)"
show("$<- on a column", q)
q <- p; q[2, "prior"] <- "normal(0, 1)"
show("[<- on a cell", q)
show("row subset", p[p$class != "cor", ])
show("rbind with a row of itself", rbind(p, p[2, ]))
show("p + set_prior()", p + set_prior("normal(0, 1)", class = "b"))
show("c(p, set_prior())", c(p, set_prior("normal(0, 1)", class = "b")))
show("as.brmsprior(p)", as.brmsprior(p))
show("as.data.frame(p)", as.data.frame(p))
show("validate_prior(sp, f)", validate_prior(set_prior("normal(0, 1)",
                                                       class = "b"), f,
                                             data = dd))
show("set_prior() + as.brmsprior(p)",
     set_prior("normal(0, 1)", class = "b") + as.brmsprior(p))
