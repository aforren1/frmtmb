# Lane meme: brms 2.23.0 on me(x) written without sdx, and frmtmb on the
# same call. Output: dev/meme-out/brms-nosdx.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib", "C:/Users/adf44/source/r/rellib-r15", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
set.seed(11)
d <- data.frame(x = rnorm(30), sx = 0.3, y = rnorm(30))
d2 <- d; d2$sx[1] <- NA
r <- function(e) tryCatch(e, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("brms standata me(x):\n"); print(r(names({attach(list(me = brms::me), name = "bme"); on.exit(detach("bme")); brms::standata(y ~ me(x), data = d)})))
cat("brms stancode me(x):\n"); print(r(substr({attach(list(me = brms::me), name = "bme2"); on.exit(detach("bme2")); brms::stancode(y ~ me(x), data = d)}, 1, 80)))
cat("frmtmb me(x):\n"); print(r(frmtmb::frm(y ~ me(x), data = d)))
cat("frmtmb NA in sdx, rows used:\n"); print(r(nobs(suppressWarnings(frmtmb::frm(y ~ me(x, sx), data = d2)))))
