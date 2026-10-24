# String and object spellings of the family, on a build.
#   Rscript dev/review-mvnlatent-strings.R <base|lane> <attach latent: yes|no>
a <- commandArgs(trailingOnly = TRUE)
base <- c("C:/Users/adf44/source/r/rellib-r11",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (a[1] == "base") base else
  c("C:/Users/adf44/source/r/wt-mvnlatent-lib", base))
suppressMessages(library(frmtmb))
if (a[2] == "yes") suppressMessages(library(frmtmb.latent))
cat("build", a[1], format(packageVersion("frmtmb")), "latent attached",
    a[2], "\n")
show <- function(label, expr) {
  r <- tryCatch({
    v <- force(expr)
    paste("OK", class(v)[1])
  }, error = function(e) paste("ERROR:", conditionMessage(e)),
  warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(sprintf("%-42s %s\n", label, substr(gsub("\n", " ", r), 1, 160)))
}
set.seed(1)
d <- data.frame(row = 1:60)
d$Y <- rbind(matrix(rnorm(60), 30), matrix(rnorm(60, 4), 30))
show('family = "mixture_mvn"', frm(bf(Y ~ 1), family = "mixture_mvn",
                                   data = d))
show('family = "Mixture_MVN"', frm(bf(Y ~ 1), family = "Mixture_MVN",
                                   data = d))
show('brmsfamily("mixture_mvn")', brmsfamily("mixture_mvn"))
show('brmsfamily("mixture_mvn", K = 2, D = 2)',
     brmsfamily("mixture_mvn", K = 2, D = 2))
show("family = mixture_mvn (the function)",
     frm(bf(Y ~ 1), family = get0("mixture_mvn",
                                  ifnotfound = function() stop("absent")),
         data = d))
show("mixture_mvn(K = 2, D = 2)", mixture_mvn(K = 2, D = 2))
show("frmtmb::mixture_mvn", frmtmb::mixture_mvn)
if (a[1] == "base") {
  us <- getFromNamespace("us_sigma", "frmtmb")
  show("base us_sigma(numeric(5), 3)", us(numeric(5), 3))
  show("base us_sigma(-400, 1, floored = TRUE)",
       print(as.matrix(us(-400, 1, floored = TRUE))))
}
