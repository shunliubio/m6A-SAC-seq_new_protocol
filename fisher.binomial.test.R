#!/usr/bin/env Rscript

start <- Sys.time()
cat("Starting Time is ",format(start,format="%Y-%m-%d %H:%M:%S"),"\n\n",sep="")

argv <- commandArgs(TRUE)

in_file <- argv[1]
bg_file <- argv[2]
out_file <- argv[3]

suppressPackageStartupMessages(library(dplyr))

in_table <- read.delim(in_file,header=T)
bg_table <- read.delim(bg_file,header=T)
in_bg_table <- left_join(in_table,bg_table[,c("motif","bg","coef")],by=c("motif"="motif"))
in_bg_table$pvalue_f <- with(in_bg_table,apply(cbind(treated_depth-treated_A_count,treated_A_count,control_depth-control_A_count,control_A_count),1,function(a) {fisher.test(matrix(a,nrow=2),alternative="greater")$p.value}))
in_bg_table$fdr_f <- p.adjust(in_bg_table$pvalue,"BH")
in_bg_table$pvalue_f <- signif(in_bg_table$pvalue,4)
in_bg_table$fdr_f <- signif(in_bg_table$fdr,4)
in_bg_table$pvalue_b <- with(in_bg_table,apply(cbind(treated_depth-treated_A_count,treated_depth,bg/100),1,function(a) {binom.test(a[1],a[2],p=a[3],alternative="greater")$p.value}))
in_bg_table$fdr_b <- p.adjust(in_bg_table$pvalue_b,"BH")
in_bg_table$pvalue_b <- signif(in_bg_table$pvalue_b,4)
in_bg_table$fdr_b <- signif(in_bg_table$fdr_b,4)

write.table(subset(in_bg_table,select = -c(bg,coef)),out_file,sep="\t",quote=F,row.names=F)

warnings()

end <- Sys.time()
cat("Ending Time is",format(end,format="%Y-%m-%d %H:%M:%S"),"\n")
runtime <- difftime(end,start,unit="mins")
cat("Used Time is",runtime,"mins\n")
