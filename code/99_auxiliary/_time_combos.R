suppressPackageStartupMessages({ library(survival) })
# 定位：到底哪个 (队列, 基因数) 组合会异常慢
zscore <- function(m){mu<-rowMeans(m,na.rm=TRUE);sd<-apply(m,1,sd,na.rm=TRUE);sd[is.na(sd)|sd<1e-8]<-1;(m-mu)/sd}
cz<-function(sv,xb){cf<-try(concordancefit(sv,xb,reverse=TRUE,std.err=FALSE),silent=TRUE)
  if(inherits(cf,"try-error")) return(NA_real_); c0<-cf$concordance
  if(is.na(c0)) return(NA_real_); if(c0<0.5) 1-c0 else c0}

for (ch in c("GSE29621","GSE39582","GSE72970","GSE17537","GSE38832")) {
  d <- readRDS(file.path("cohorts", paste0(ch, ".rds")))
  z <- zscore(d$expr); cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), ]
  tm <- as.numeric(cl$time); st <- as.integer(cl$status); kp <- !is.na(tm) & !is.na(st)
  z <- z[, kp]; tm <- tm[kp]; st <- st[kp]
  sv <- Surv(tm, st); nk <- seq_len(nrow(z)); n <- length(tm)
  for (ng in c(16, 27, 35, 49, 50)) {
    set.seed(1); xb <- matrix(0, n, 1)
    t0 <- Sys.time()
    for (b in 1:300) { xb[,1] <- colMeans(z[sample(nk, ng), , drop=FALSE]); cz(sv, xb) }
    el <- as.numeric(Sys.time() - t0, units = "secs")
    cat(sprintf("%-11s n=%3d  300次 %6.2f 秒  -> B=3000 预计 %6.1f 分钟\n",
                ch, ng, el, el*10/60))
  }
  rm(d, z); invisible(gc(FALSE))
}
