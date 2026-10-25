## brms 2.23.0's coding of an integer mo() predictor that starts at 1:
## Jmo and the first values of Xmo, beside the table of k. Data seed 1024
## (dev/brmsnames-1024-cases.R).
##   Rscript dev/brmsnames-1024-probe1.R > dev/brmsnames-1024-log/probe1.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("dev/brmsnames-1024-cases.R")
d <- brmsnames_data()
s <- suppressMessages(brms::standata(brms::bf(y ~ mo(k)), data = d))
print(s$Jmo); print(table(d$k)); print(head(s$Xmo_1))
