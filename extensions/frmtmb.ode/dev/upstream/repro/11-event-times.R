## deSolve inserts event times that are not already in 'times', so the solution
## matrix it returns has MORE rows than 'times'.  Check what RTMBode does with
## that, both before and after the truncation guard.
##
## Usage: Rscript 11-event-times.R <lib>

args <- commandArgs(TRUE)
if (length(args) >= 1 && nzchar(args[1])) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages({ library(RTMB); library(RTMBode) })

pk <- function(t, y, p) list(c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2]))
y0 <- c(depot = 100, central = 0)
times <- seq(0, 24, by = 2)      # 13 output times

run <- function(dosetimes, label) {
    ev <- list(data = data.frame(var = "depot", time = dosetimes,
                                 value = 50, method = "add"))
    cat("\n== ", label, " (dose times ",
        paste(dosetimes, collapse = ", "), ") ==\n", sep = "")
    n <- suppressWarnings(nrow(deSolve::ode(y0, times, pk, c(1.2, 0.35), events = ev)))
    cat("length(times) =", length(times), "  nrow(numeric solve) =", n, "\n")
    f <- function(pl) sum(RTMBode::ode(y0, times, pk, exp(pl$par), events = ev)[, "central"])
    r <- suppressWarnings(tryCatch({
        obj <- MakeADFun(f, list(par = log(c(1.2, 0.35))), silent = TRUE)
        p <- log(c(1.2, 0.35))
        list(v = obj$fn(p), g = as.numeric(obj$gr(p)))
    }, error = function(e) conditionMessage(e)))
    if (is.character(r)) { cat("ERROR:", gsub("\\s+", " ", r), "\n"); return(invisible()) }
    h <- 1e-5
    fd <- sapply(1:2, function(i) {
        pp <- pm <- log(c(1.2, 0.35)); pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
        num <- function(q) {
            s <- suppressWarnings(deSolve::ode(y0, times, pk, exp(q), events = ev))
            sum(s[s[, 1] %in% times, "central"])
        }
        (num(pp) - num(pm)) / (2 * h)
    })
    cat("value :", signif(r$v, 10), "\n")
    cat("AD    :", signif(r$g, 8), "\n")
    cat("FD    :", signif(fd, 8), "\n")
    cat("max rel err:", signif(max(abs(r$g - fd) / pmax(abs(fd), 1e-8)), 3), "\n")
}

run(c(6, 12, 18), "event times already in 'times'")
run(c(5, 11, 17), "event times NOT in 'times'")
