.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
f <- function(p) 0 * log(p[1]) + p[1]^2
tp <- RTMB::MakeTape(f, 2)
print(suppressWarnings(f(-1))); print(tp(-1))
g <- function(p) { u <- p[1] + abs(p[1]); u / u - 1 + p[1]^2 }
print(RTMB::MakeTape(g, 2)(-1))
