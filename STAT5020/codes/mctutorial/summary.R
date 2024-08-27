load("bflinear-20230322-16.RData")
bf21 = resultlst

load("bfalter-20230322-16.RData")
bf13 = resultlst

load("dictrue-20230322-15.RData")
dic1 = resultlst

load("diclinear-20230322-15.RData")
dic2 = resultlst

load("dicalter-20230322-15.RData")
dic3 = resultlst


{
  for (i in 1:10) {
    print(paste(i, "&", bf21$lbf[i], "&", bf13$lbf[i], "&", dic1$dic[i], "&", dic2$dic[i], "&", dic3$dic[i], "\\"))
  }
  print(paste("Mean", "&", bf21$lbf.mean, "&", bf13$lbf.mean, "&",dic1$dic.mean, "&", dic2$dic.mean, "&", dic3$dic.mean, "\\"))
  print(paste("SD", "&", bf21$lbf.mean, "&", bf13$lbf.sd, "&",dic1$dic.sd, "&", dic2$dic.sd, "&", dic3$dic.sd, "\\"))
}
