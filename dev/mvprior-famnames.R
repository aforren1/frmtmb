# Lane wt-mvprior: the family name and location dpars of each family
# whose location is several dpars, to key the location rule on.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressMessages(library(frmtmb))
show <- function(f) {
  cat(sprintf("%-28s family=%-18s primary=%s\n", deparse(substitute(f)),
              f[["family"]], paste(f[["primary_dpars"]], collapse = ",")))
}
show(categorical())
show(multinomial(K = 3))
if (exists("dirichlet")) show(dirichlet())
if (exists("logistic_normal")) show(logistic_normal())
show(mixture(gaussian(), gaussian()))
show(mixture(gaussian(), poisson()))
show(mixture_mvn(K = 2, D = 2))
cat("mix slot on mixture():", !is.null(mixture(gaussian(), gaussian())[["mix"]]),
    " on mixture_mvn():", !is.null(mixture_mvn(K = 2, D = 2)[["mix"]]), "\n")
