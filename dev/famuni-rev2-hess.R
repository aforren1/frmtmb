.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
g <- frmtmb:::nan_unless_pos
tp <- RTMB::MakeTape(function(p) g(p), 1)
hp <- RTMB::MakeTape(function(p) tp$jacfun()(p), 1)
for (e in c(-100, -140, -150, -153, -154, -155, -160, -200, -300, -307, -308)) {
  x <- 10^e
  cat(sprintf("x = 1e%d: grad %s hess %s\n", e, format(tp$jacobian(x)),
              format(hp$jacobian(x))))
}
