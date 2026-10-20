# Mentions of a lane name, the word "lane", or a `dev/` path in the
# 0.72.0 section of each package's NEWS (dev/rel072-news-strip.R's
# mentions()).
#   Rscript dev/rel072-news-count.R
source("dev/rel072-news-strip.R")
for (f in c("NEWS.md", "extensions/frmtmb.sample/NEWS.md",
            "extensions/frmtmb.eam/NEWS.md", "extensions/frmtmb.latent/NEWS.md",
            "extensions/frmtmb.coupling/NEWS.md")) {
  x <- sub("\r$", "", readLines(f, warn = FALSE))
  h <- grep("^# ", x)
  cat(sprintf("%-36s %s mentions: %d\n", f, x[h[1]],
              mentions(x[h[1]:(h[2] - 1L)])))
}
