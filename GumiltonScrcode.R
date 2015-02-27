
### Plot MA for multiple factor designed Fit
MultiPlotMA <- function(fit, pcut=0.05){
  for (i in 1:ncol(fit)){
    png(paste("MA Plot",colnames(fit)[i],".png"), width = 700, height = 650, units = "px")
    dat <- topTable(fit,coef=i, number=Inf,adjust.method="BH",sort.by = "P")
    plot(dat$AveExpr, dat$logFC, main=paste("MA Plot",colnames(fit)[i]),
         xlab="Average Log-Expression", ylab="Log2 Fold Change",
         col=ifelse(dat$adj.P.Val< pcut,"red","black"), pch=19)
    abline(h=c(-1,0,1),col=c("yellow","black","yellow"))
    dev.off() 
  } 
}

### Plot MA for multiple factor with contrast designed Fit2

multi_plotMA <- function(fit, p =.1){
  contrasts <- fit$contrasts
  vs <- c()
  for (i in 1:ncol(contrasts)){
    x <- contrasts[,i]
    vs <- c(vs, paste(paste(names(x)[x>0],collapse=""),"vs",paste(names(x)[x<0],collapse="")))
  }
  for (i in 1:ncol(contrasts)){
    dat <- topTable(fit,coef=i, number=Inf, adjust.method="BH", sort.by = "none")
    png(paste("MA Plot",vs[i],".png"), width = 700, height = 650, units = "px")
    plot(dat$AveExpr, dat$logFC, main=paste("MA Plot",vs[i]), pch=19, cex=.7,
         col = ifelse(dat$adj.P.Val< p, "red","black"), xlab="Average Expression",
         ylab = "Log2 Fold Change")
    abline(h=c(-1,0,1),col=c("yellow","black","yellow"))
    dev.off()
  }
}


### Density Plot

DensityPlot <- function(data.frame,main=""){
  color <- rainbow(1024)
  xlim_vec <- range(density(data.frame)$x)
  
  y_density_range <- function(data){
    y<-range(density(data)$y)
    return (c(y))
  }  
  
  ylim_vec <-range(apply(data.frame,2,y_density_range))
  
  plot.new()
  plot.window(xlim_vec, ylim_vec)
  
  apply(data.frame,2,
        function(x){points(density(x),
                           type = "l", lty=3,
                           ylim = ylim_vec,
                           xlim = xlim_vec,
                           main="",col=sample(color,1))}
  )
  box()
  axis(1)
  axis(2)
  title(main)
}

### Bar plot for genes in a genelist

genesbarplot <- function(genelist, data){  
  for(i in 1:length(genelist)){
    png(filename=paste(genelist[i],"Expressions.png"), 
        width=10, height=10,units="in", res=100)
    par(mar = c(15,5,5,5))
    barplot(as.matrix(data[genelist[i],]),las=2,cex.names = 1.2)
    title(paste(genelist[i],"Expressions"), cex.main = 2)
    dev.off()
  }  
}

### Plot survival figures for individual gene

plotSurv <- function(BCR_dat, BCR_info, geneset, dir="./", plot=T){
  library(xlsx)
  library(plyr)
  dir.create(file.path(dir), showWarnings = FALSE,recursive =T)
  sig_gene <- data.frame(gene=c(),pval=c())
  for(i in geneset){
    if(!i %in% rownames(BCR_dat)) next
    
    temp <- t(BCR_dat[i,])
    gene <- colnames(temp)
    temp <- as.data.frame(temp)
    temp$Sample.ID <- rownames(temp)
    temp <- mutate(temp,label = ifelse(temp[,1] < median(temp[,1]),
                                       "Low","High"))
    dat_surv <- merge(BCR_info,temp, by = "Sample.ID")
    
    fit <- survfit(Surv(BCR_FreeTime, BCR_Event) ~ label, data = dat_surv) 
    dif <- survdiff(Surv(BCR_FreeTime, BCR_Event) ~ label, data = dat_surv)
    pval <- round(1 - pchisq(dif$chisq, length(dif$n) - 1),3)
    
    if (plot){
      png(paste0(dir,"/Taylor2010CancerCell_MKSCC_",gene, 
                 "_BCRFree_Survival.png", sep=""),
          width=13, height=10,units="in", res=100)
      par(mar = c(5,5,5,5), cex=2)
      
      plot(fit, lty = 1:2, xlab = "BCR FreeTime", ylab = "BCR Free Surviving", 
           lwd = 2, col=c("red","blue"))
      
      legend(x = 10, y = 0.28, legend = c(paste0(gene," High ",dif$obs[1],"/",dif$n[1]),
                                          paste0(gene," Low ",dif$obs[2],"/",dif$n[2])),
             lty = c(1:2), lwd = 2, cex = 0.8, col=c("red","blue"))
      #text(0.05, 0.05, "(A) TCGA", pos = 4)
      text(10, 0.4, paste("P =",pval), pos = 4)
      title(paste("BCR Free Survival on",gene))
      dev.off()
      
    }
    sig_gene <- rbind(sig_gene,cbind(gene, pval))
  }
  sig_gene$adjp <- p.adjust(sig_gene$pval, method="BH")
  
  return((sig_gene))
  write.table(sig_gene, paste0(dir,"/Summary_P_Value_Taylor2010CancerCell_MKSCC_",
                               substitute(geneset), "_BCRFree_Survival.txt"), row.names=F)
}

### Read expression raw data for multiple samples

readfile <- function(path){
  li <- grep("Genes_ReadCount.txt",list.files(directory),value=TRUE)
  dat <- list()
  for(i in 1:length(li)){
    
    print(paste("reading",li[i]))
    dat <- c(dat,read.delim(paste0(directory,li[i]),header=F))
    names(dat)[2*i] <- gsub("Sample_(.*)_Genes_ReadCount.txt","\\1",li[i])
  }
  dat <- as.data.frame(dat)
  rownames(dat) <- dat$V1
  dat <- dat[,-seq(1,2*i,2)]
  return(dat)
}

### Separate list of genes into list of two with up and down by the foldchange

sep_up_down <- function(list, data){
  down <- list[data[list,"logFC"]<0]
  up <- list[data[list,"logFC"]>0]
  return(list(up=up,down=down))  
}

### Convert tab-deliminated file .txt to .xlsx

txt2xlsx <- function(path){
  library(xlsx)
  li <- grep(".txt",list.files(path),value=TRUE)
  xlsx <- grep(".xlsx",list.files(path),value=TRUE)
  print(li)
  for (i in li){
    if(! substr(i,1,nchar(i)-3) %in% substr(xlsx,1,nchar(xlsx)-4)){
      write.xlsx(read.delim(paste0(path,i)), paste0(path,substr(i,1,nchar(i)-3), "xlsx"))
    }
  }
}