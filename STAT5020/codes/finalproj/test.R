library(RColorBrewer)
library(colorjam)
library(npreg)

work.path = "D:/OneDrive/OneDrive - sjtu.edu.cn/Records/CUHK/Course/STAT5020/codes/finalproj"
rds.path = "D:/Project/Database/ABCD2.0/nda2.0.1_orig.Rds"
abcd_repo = readRDS(rds.path)
eventname = abcd_repo$eventname
unique(eventname)
abcd_repo.baseline = abcd_repo[abcd_repo$eventname == "baseline_year_1_arm_1",]

saveRDS(abcd_repo.baseline, file = "D:/Project/Database/ABCD2.0/nda2.0.1_baseline.Rds")
save(abcd_repo.baseline, file = "D:/Project/Database/ABCD2.0/nda2.0.1_baseline.RData")

var.names = colnames(abcd_repo)
smp.names = abcd_repo.baseline$src_subject_id
"cbcl_scr_dsm5_opposit_t" %in% var.names
length(unique(smp.names))
"sex" %in% var.names
"ethnity" %in% var.names
write.table(var.names, file = "var_names.txt")
"sleepdisturb1_p" %in% var.names
sleep.vars = read.table("sleep.txt", header = F, stringsAsFactors = FALSE)$V1
sleep.df = abcd_repo.baseline[,sleep.vars]
sleep.df.n = sapply(sleep.df, as.numeric)
sleep.cor = cor(sleep.df.n, method = "spearman")
sleep.sd = apply(sleep.df.n,2, sd)

colmap = brewer.pal(n = 7, name ="Greens")
color.name = closestRcolor(number2color(sleep.sd, colmap, ncol = 20))

heatmap(sleep.cor, RowSideColors = color.name)


# selecting col
library(readxl)
col.sel = read_excel("selectedcols.xlsx")
abcd_sel.baseline = abcd_repo.baseline[,col.sel$var.name]
rownames(abcd_sel.baseline) = abcd_repo.baseline$src_subject_id
saveRDS(abcd_sel.baseline, file = "abcd_sel.baseline.Rds")

rsfqc = abcd_repo.baseline$fsqc_qc == "1" & as.numeric(abcd_repo.baseline$rsfmri_cor_network.gordon_ntpoints) > 375
abcd_sel.logical = is.na(abcd_sel.baseline) | abcd_sel.baseline == ""
abcd_sel.cleaned = abcd_sel.baseline[rowSums(abcd_sel.logical) == 0 & rsfqc,]
saveRDS(abcd_sel.cleaned, file = "abcd_sel.cleaned.Rds")

abcd_sel.cleaned.n = as.data.frame(apply(abcd_sel.cleaned, 2, as.numeric))
abcd_sel.cleaned.n[,"sex"] = abcd_sel.cleaned$sex
rownames(abcd_sel.cleaned.n) = rownames(abcd_sel.cleaned)
colnames(abcd_sel.cleaned.n) = col.sel$y

# convert dataframe to WinBUGS format
N = dim(abcd_sel.cleaned.n)[1]
ohc1 = as.data.frame(matrix(rep(0, 5 * N), nrow = N, ncol = 5))
for(i in 1:N) {
  ohc1[i,abcd_sel.cleaned.n$c1[i]] = 1
}
colnames(ohc1) = paste0("c1.", 1:5)
abcd_sel.cleaned.onehot = cbind(abcd_sel.cleaned.n[,1], 
                                abcd_sel.cleaned.n[,3], 
                                ohc1[,2:5], 
                                abcd_sel.cleaned.n[,4:9],
                                abcd_sel.cleaned.n[,16:dim(abcd_sel.cleaned.n)[2]],
                                abcd_sel.cleaned.n[,10:15] # change location
                                )
colnames(abcd_sel.cleaned.onehot) = c("g", paste0("d",1:5), paste0("z",1:17))

for (i in 1:17) {
  hist(abcd_sel.cleaned.onehot[,paste0("z",i)], main = paste0("z",i))
}

abcd_sel.cleaned.onehot[,paste0("z",1:6)] = sapply(
  abcd_sel.cleaned.onehot[,paste0("z",1:6)], 
  function(v){return(as.integer(v <= 54)*1 + (v >= 55 & v <= 64)*2 + (v >= 65)*3)}
)

abcd_sel.cleaned.onehot[,paste0("z",7:11)] = sapply(
  abcd_sel.cleaned.onehot[,paste0("z",7:11)],
  as.integer
  )

abcd_sel.cleaned.onehot$d1 = abcd_sel.cleaned.onehot$d1 - 120

for (i in 1:17) {
  hist(abcd_sel.cleaned.onehot[,paste0("z",i)], main = paste0("z",i))
}

sum(abcd_sel.cleaned.onehot$g=="F")
sum(abcd_sel.cleaned.onehot$g=="M")

saveRDS(abcd_sel.cleaned.onehot, file = "abcd_sel.cleaned.onehot.Rds")
write.csv(abcd_sel.cleaned.onehot, file = "abcd_sel.cleaned.onehot.csv")


