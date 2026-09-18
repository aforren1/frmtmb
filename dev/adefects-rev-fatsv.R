setwd("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-log")
f <- function(a) read.delim(sprintf("falsealarm-%s.tsv", a), quote = "",
                            stringsAsFactors = FALSE)
b <- f("base"); l <- f("lane")
cat("base rows", nrow(b), " lane rows", nrow(l), "\n")
for (k in unique(c(b$kind, l$kind))) {
  for (nm in c("base", "lane")) {
    x <- if (nm == "base") b else l
    x <- x[x$kind == k, ]
    u <- x[!duplicated(x$text), ]
    cat(sprintf("%-8s %-4s rows=%3d distinct=%3d both_accept=%3d both_refuse=%3d frmRef_brmsAcc=%3d frmAcc_brmsRef=%3d\n",
                k, nm, nrow(x), nrow(u), sum(u$brms & u$frm),
                sum(!u$brms & !u$frm), sum(u$brms & !u$frm),
                sum(!u$brms & u$frm)))
  }
}
cmp <- function(k) {
  xb <- b[b$kind == k, ]; xl <- l[l$kind == k, ]
  xb <- xb[!duplicated(xb$text), ]; xl <- xl[!duplicated(xl$text), ]
  m <- merge(xb, xl, by = "text", suffixes = c(".b", ".l"))
  ch <- m$frm.b != m$frm.l
  cat(sprintf("%s merged=%d changed=%d  [%s]\n", k, nrow(m), sum(ch),
              paste(m$text[ch], collapse = " | ")))
  cat(sprintf("  base-only=%d lane-only=%d\n",
              sum(!xb$text %in% xl$text), sum(!xl$text %in% xb$text)))
}
for (k in unique(c(b$kind, l$kind))) cmp(k)
