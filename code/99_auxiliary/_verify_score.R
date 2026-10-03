suppressPackageStartupMessages({library(survival)})
# 用最朴素的方式实现同一件事，逐个和向量化版本对比
zscore <- function(m){mu<-rowMeans(m,na.rm=TRUE);sd<-apply(m,1,sd,na.rm=TRUE);sd[is.na(sd)|sd<1e-8]<-1;(m-mu)/sd}
eval(parse(file="/tmp/cox_score_fn.R"))

d <- readRDS("cohorts/GSE17536.rds")
z <- zscore(d$expr); cl <- d$clin; cl <- cl[match(colnames(z), cl$sample), ]
tm <- as.numeric(cl$time); st <- as.integer(cl$status); kp <- !is.na(tm) & !is.na(st)
z <- z[, kp]; tm <- tm[kp]; st <- st[kp]
cat("样本", length(tm), "事件", sum(st), "\n")
cat("时间重复情况：", length(unique(tm)), "个不同时间 /", length(tm), "\n")

# -------- 朴素参考实现（直接照定义循环） --------
brute_score <- function(x, time, status) {
  et <- unique(time[status == 1])
  U <- 0; V <- 0
  for (tt in et) {
    R  <- which(time >= tt)            # 风险集
    nb <- sum(time == tt & status == 1) # 该时点事件数
    xb <- mean(x[R])
    for (k in which(time == tt & status == 1)) U <- U + (x[k] - xb)
    # Breslow 信息
    V  <- V + nb * (mean(x[R]^2) - xb^2)
  }
  c(U = U, V = V, z = if (V > 0) U / sqrt(V) else 0)
}

set.seed(11)
gs <- sample(rownames(z), 6)
zz <- cox_score(z, tm, st)
tb <- t(vapply(gs, function(g) brute_score(as.numeric(z[g, ]), tm, st), numeric(3)))
vv <- data.frame(gene = gs,
                 U_brute = round(tb[, "U"], 4), U_vec = NA,
                 V_brute = round(tb[, "V"], 4), V_vec = NA,
                 z_brute = round(tb[, "z"], 4), z_vec = round(zz[gs], 4),
                 check.names = FALSE)

# 从向量化版本内部取 U / V（重算一次，把中间量打印出来）
debug_cox <- function(z, time, status) {
  ord <- order(time, decreasing = TRUE)
  X <- z[, ord, drop = FALSE]; t_o <- time[ord]; s_o <- status[ord]
  r <- rle(t_o); ends <- cumsum(r$lengths)
  ne <- vapply(seq_along(ends), function(k) {
    ii <- if (k == 1L) 1:ends[1] else (ends[k-1] + 1):ends[k]
    sum(s_o[ii] == 1)
  }, numeric(1))
  keep <- ne > 0
  ends <- ends[keep]; ne <- ne[keep]
  CX  <- t(apply(X, 1, cumsum))
  CX2 <- t(apply(X * X, 1, cumsum))
  SXblk <- t(apply(X * (s_o == 1), 1, cumsum))
  U <- numeric(nrow(z)); V <- numeric(nrow(z))
  for (i in seq_along(ends)) {
    eu <- ends[i]
    S1 <- CX[, eu] / eu
    Se <- SXblk[, eu] - if (i == 1L) 0 else SXblk[, ends[i-1]]
    S2 <- CX2[, eu] / eu
    U <- U + Se - ne[i] * S1
    V <- V + ne[i] * pmax(S2 - S1 * S1, 1e-8)
  }
  list(U = setNames(U, rownames(z)), V = setNames(V, rownames(z)))
}
dbg <- debug_cox(z, tm, st)
vv$U_vec <- round(dbg$U[gs], 4)
vv$V_vec <- round(dbg$V[gs], 4)
print(vv)
cat("\nU 相关性:", round(cor(tb[, "U"], dbg$U[gs]), 4),
    "  V 相关性:", round(cor(tb[, "V"], dbg$V[gs]), 4), "\n")
