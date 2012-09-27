rm(list=ls())

library(plyr)
library(reshape2)
library(abind)
library(xtable)


load.file <- "~/Documents/R_packages/docs/trustOptim/test1.Rdata"
save.file <- "~/Documents/R_packages/docs/trustOptim/compare_summary.Rdata"
tex.file <- "~/Documents/R_packages/docs/trustOptim/compare_summary.tex"

load(load.file)

x <- melt(res)
x$stat <- reorder(x$stat,new.order=c("time","nrm.gr","max.abs.gr","iters"))
x$method <- reorder(x$method,new.order=c("SparseFD","Sparse","BFGS.optim","CG","trust"))


y1 <- acast(x[x$k==2,], method~stat~N+k,fun.aggregate=mean,value.var="value")
y2 <- acast(x[x$k==10,], method~stat~N+k,fun.aggregate=mean,value.var="value")
y <- abind(y1,y2,along=2)
z <- NULL
for (i in 1:3) {
##  z <- rbind(z,cbind(rownames(y[,,i]),i,y[,,i]))
    z <- rbind(z,cbind(i,y[,,i]))
}

save(x,y,z,file=save.file)

tab <- xtable(z,
              display=c("s","d","f","g","g","f","f","g","g","f"),
              digits=c(1,3,3,3,3,3,3,3,3,3),
              rownames=TRUE
              )
print(tab,file=tex.file)



