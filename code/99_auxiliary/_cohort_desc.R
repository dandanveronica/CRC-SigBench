d <- readRDS("cohorts/TCGA_COADREAD.rds")
cat("TCGA: n=", nrow(d$clin), " genes=", nrow(d$expr), " events=", sum(d$clin$status==1,na.rm=TRUE), "\n")
cat(" time range:", range(d$clin$time,na.rm=TRUE), "月\n")
