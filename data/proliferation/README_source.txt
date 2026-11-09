# 增殖基因集来源留痕（post-hoc 分析 ⑦ 用）
# 
# 基因集：KEGG hsa04110 "Cell cycle"，124 个人类基因符号
# 获取：Enrichr REST API，library=KEGG_2021_Human，term="Cell cycle"
#       （2026-11-08 抓取；KEGG 官网被本机网络拦截，改用 Enrichr 镜像，
#        该库逐字对应 KEGG hsa04110 pathway）
# 引用建议：Kanehisa & Goto 2000 (KEGG)；稿件中写 "the KEGG cell-cycle pathway
#           (hsa04110, n=124 genes)"
# 用途：逐患者增殖评分 = 队列内 124 基因表达 z-score 的等权平均（与签名打分同口径）
