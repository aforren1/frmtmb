source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
nss_report_env()

cat("\nStanHeaders:", format(packageVersion("StanHeaders")), "\n")

nm <- getNamespaceExports("RTMB")
cat("\nRTMB exports matching value/context/tape:\n")
print(sort(grep("value|Value|context|Tape|tape", nm, value = TRUE)))

cat("\nRTMB internal names matching value/context:\n")
inm <- ls(asNamespace("RTMB"), all.names = TRUE)
print(sort(grep("value|Value|ad_context|getVal", inm, value = TRUE)))
