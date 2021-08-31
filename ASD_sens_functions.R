flattenAnn <- function(array.type)
  # flatten 450k or EPIC array annotation
  # Belinda Phipson
  # 10 February 2016
  # Updated 7 July 2016
{
  library(stringr)
  library(limma)
  
  if(array.type=="450K")    
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  else
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
  # get rid of the non-CpG sites
  strlen<-str_length(rownames(anno))
  ann.keep<-anno[strlen==10,]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  
  flat$cpg<- substr(rownames(flat),1,10)
  
  flat$alias <- alias2SymbolTable(flat$symbol)
  
  eg <- toTable(org.Hs.egSYMBOL2EG)
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  
  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
}




betavsbeta <- function(res1,res2,title,folder){
  pdf(folder)
    x <- res1[,1]*100
    y <- res2[rownames(res1),1]*100
    smoothScatter(x,y, xlim=c(-30,30), ylim=c(-30,30),
                  bandwidth=0.01,colramp=colorRampPalette(c('white','slateblue3','slateblue4')),
                  main=title,xlab="SV",ylab="Known Covariates")
    abline(h=0);abline(v=0)
    legend("bottomleft",paste0("r=",round(cor(x,y),2)),cex=3,bty='n')
  dev.off()
  
}

pvp <- function(res1,res2,title,folder){
  pdf(folder)
  x <- -log(res1[,'P.Value'],10)
  y <- -log(res2[rownames(res1),'P.Value'],10)
  smoothScatter(x,y, xlim=c(0,10), ylim=c(0,10),
                bandwidth=0.01,colramp=colorRampPalette(c('white','slateblue3','slateblue4')),
                main=title,xlab="SV",ylab="Known Covariates")
  abline(h=0);abline(v=0)
  legend("bottomleft",paste0("r=",round(cor(x,y),2)),cex=3,bty='n')
  dev.off()
  
}



###########################################################
### formats top cpg hits table
###########################################################
#
# sshits = output from toptable after using lmFit from limma
# what = vector of annotation info to include:
# array = 450k or EPIC
#

fmt.tophits <- function(sshits, array=c('450k','EPIC')){
  #check if sshits is what expected
  if(!identical(names(sshits),c('logFC','AveExpr','t','P.Value','adj.P.Val','B'))){
    stop('Columns do not match the expected format from topTable() output.')
  }

  #default assumption 450k
  array <- array[1]
  
  #annotation info
  sshits <- append.annotate(sshits, what=c('chr','pos','gene'), array=array)
  
  #standard error
  sshits$SE <- sshits$logFC/sshits$t
  
  #keep effect, p-val, fdr, chr, pos, gene columns and format
  sshits <- sshits[,c(9,7,8,2,1,4,10,5)]
  names(sshits) <- c('GeneSymbol','Chr','Pos','MeanMeth','DiffMeth','P','SE','FDR')
  sshits$MeanMeth <- round(sshits$MeanMeth*100,1)
  sshits$DiffMeth <- round(sshits$DiffMeth*100,2)
  sshits$SE <- round(sshits$SE*100,2)
  sshits$P <- signif(sshits$P,digits=2)
  sshits$GeneSymbol <- ifelse(sshits$GeneSymbol=='','-',sshits$GeneSymbol)
  return(sshits)
}



###########################################################
### appends cpg annotation columns to a data frame
###########################################################
#
# X = data frame with cpg names as row names or in a specified column
# what = vector of annotation info to include:
#   chr - chromosome
#   pos - position
#   gene - USCS annotated gene
#   island - relation to islands
# array = 450k or EPIC
# cpg = column name for cpgs, if its not the rownames of X
#
append.annotate <- function(X, what=c('chr','pos','gene','island'), array='450k', cpg=NULL){
  
  #load appropriate annotation package
  if(array=='450k' | array=='450K'){
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }else if(array=='epic' | array=='EPIC'){
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  }else{
    stop('Invalid array specified. 450k or EPIC only.')
  }
  
  data(Locations)
  data(Other)
  data(Islands.UCSC)
  
  #check X rownames
  data(Locations)
  rownamen <- rownames(X) #save original rownames just in case they are actually important
  if(is.null(cpg)){
    check <- rownames(X)[rownames(X) %in% rownames(Locations)]
  }else{
    if(cpg %in% colnames(X)){
      check <- X[X[,cpg] %in% rownames(Locations),cpg]
      rownames(X) <- X[,cpg]
    }else{
      stop('The column name you gave for CpG names is not in the dataframe.')
    }
  }
  
  if(length(check) < nrow(X)){
    stop(paste0(nrow(X)-length(check), ' of the CpG names are not in the annotation file, something is wrong, check your rownames or cpg column.'))
  }
  
  #trim and order the annotaiton files
  Locations <- Locations[rownames(X),]
  Other <- Other[rownames(X),]
  Islands.UCSC <- Islands.UCSC[rownames(X),]
  
  #add things to data frame
  if('chr' %in% what){
    X$chr <- Locations$chr
  }
  if('pos' %in% what){
    X$pos <- Locations$pos
  }
  if('gene' %in% what){
    X$gene <- Other$UCSC_RefGene_Name
  }
  if('island' %in% what){
    X$island <- Islands.UCSC$Relation_to_Island
  }
  
  rownames(X) <- rownamen
  return(X)
}



###############################################################################################
# makes sva variables
# beta - beta matrix
# pd - phenotype file
# var - variable of interest
#
# returns your pd back, with sv variables merged
###############################################################################################
sva.fit <- function(beta,pd,var){
  nullmod <- model.matrix(~1,data=pd)
  fullmod <- model.matrix(~pd[,var])
  sva.fit <- sva(beta, fullmod, nullmod)
  sv <- data.frame(sva.fit$sv)
  colnames(sv) <- paste('sv',1:ncol(sv),sep='')
  pd.sv <- cbind(pd,sv)
  
  return(pd.sv)
}

###############################################################################################
# heatmap of sv association with variables
# pd - phenotype file, with sv vars with names sv1, sv2, etc.
# covar - vector of known variables of interest
#
###############################################################################################
bring.the.heatmap.sv <- function(pd,covar,folder){
  library(gplots)
  
  nsv <- length(grep('sv',names(pd)))
  
  #generate correlation matrix between SVs and variables
  correlations <- matrix(0,nrow=length(covar),ncol=nsv)
  colnames(correlations) <- paste('sv',1:nsv,sep='')
  rownames(correlations) <- covar
  
  #take svs
  svs <- pd[,colnames(correlations)]
  
  #correlations between svs and variables:
  # cor test if continuous
  # anova is multi category
  # t-test if two category
  for(col in covar){
    var <- pd[,col]
    #check to see if it is numeric and isn't really a two category thing
    if(class(var)=='numeric' & dim(table(var))>2){
      correlations[col,] <- apply(svs,2,FUN = function(x,y) cor.test(x,y)$p.value,y=var)
    }else{
      if(length(levels(factor(var)))>2){
        correlations[col,] <- apply(svs,2,FUN = function(x,y) anova(lm(x~y))$'Pr(>F)'[1],y=var)
      }else{
        correlations[col,] <- apply(svs,2,FUN = function(x,y) t.test(x~y)$p.value,y=var)
      }
    }
  }
  
  pdf(paste0(folder,'heatmap-sv.pdf'))
  color.gradient <- colorRampPalette(c('firebrick','ivory3','royalblue'),bias=3)
  heatmap.2(correlations,Rowv=FALSE,Colv=FALSE,cexRow=1,cexCol=1,dendrogram='none',density.info='none',trace='none',
            col=color.gradient(20), key.xlab='p-value', main='Associations Between SVs\n and Given Covariates',
            lmat=rbind(c(0,3,3,3),c(2,0,4,0),c(0,1,1,1)), lhei=c(0.5,1,3.5), lwid=c(0.12,0.145,0.4,0.335), margin=c(4,7))
  dev.off()
  
  return('check your plot')
}

###############################################################################################
# makes plots to see effect on lambda from sequential addition
# beta - beta matrix
# pd - phenotype file with svs
# var - variable of interest
# folder - filepath where you want stuff saved
###############################################################################################
sva.sequential <- function(beta,pd,var,folder){
  sv.vars <- as.formula(paste0("~pd$",var))
  nsv <- length(grep('sv',names(pd)))
  
  #set up object to save lambda values
  lambdas <- matrix(0,nrow=1,ncol=nsv+1)
  colnames(lambdas) <- paste('use',0:nsv,'sv',sep='')
  
  #set up data frame to save beta coefficients
  mega.table <- data.frame(matrix(0,nrow=nrow(beta),ncol=nsv+1))
  colnames(mega.table) <- paste('use',0:nsv,'sv',sep='')
  rownames(mega.table) <- rownames(beta)
  
  for(numsvs in 0:nsv){
    design <- model.matrix(sv.vars)
    fit <- lmFit(beta,design)
    fit <- eBayes(fit)
    res.sv <- topTable(fit,coef=2,nrow(beta))
    
    #add in lambdas and betas
    lambdas[numsvs+1] <- lambda(res.sv)
    mega.table[,numsvs+1] <- res.sv[rownames(beta),'logFC']
    
    png(paste0(folder,'QQ-',numsvs,'SVs.png'))
      qq(res.sv$P.Value,  xlab='Expected -log10(p)', 
            ylab='Observed -log10(p)',main=paste(numsvs,'SVs','lambda =',round(lambdas[numsvs+1],2)))
    dev.off()

    sv.vars <- as.formula(paste('~',as.character(sv.vars)[2],' + ', 'pd$sv', (numsvs+1),sep=''))
  }
  
  #lambda vs number of svs
  pdf(paste0(folder,'lambdas-surrogate-vars.pdf'))
    plot(0:nsv,lambdas[1,],xlab='Num SVs',ylab='Lambda');abline(h=1)
  dev.off()
  
  #correlation plot of beta vs beta
  M <- cor(mega.table)
  pdf(paste0(folder,'corrplot-sv-betas.pdf'))
    corrplot(M,type='upper',method='ellipse',addCoef.col='black',diag=F)
  dev.off()
  
  save(mega.table,file=paste0(folder,'mega.table.rda'))
  
  return(lambdas)
}


###############################################################################################
# input the full toptable() from limma results, get a lambda inflation factor back
###############################################################################################
lambda <- function(ss.hits){
  qchisq(median(ss.hits$P.Value,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
}


###############################################################################################
# restrict to specific genomic region in relation to cpg islands
###############################################################################################
genomic.region <- function(X, region, anno='450k'){
  #X = anything with cpgs as rownames
  #region = what area you want
  #anno = 450k or epic
  
  if(anno=='450k' | anno=='450K'){
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }else if(anno=='epic' | anno=='EPIC'){
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  }else{
    stop('Dunno that array m8, 450k or EPIC only.')
  }
  
  #pick out CpGs in X that are in region
  data(Islands.UCSC)
  if(region=='Shore' | region=='Shelf'){
    #if only 'Shore' or 'Shelf' are given, north and south are combined
    Islands.UCSC <- Islands.UCSC[rownames(X),]
    Islands.UCSC <- Islands.UCSC[Islands.UCSC$Relation_to_Island==paste('N_',region,sep='') | Islands.UCSC$Relation_to_Island==paste('S_',region,sep=''),]
  }else{
    Islands.UCSC <- Islands.UCSC[rownames(X),]
    Islands.UCSC <- Islands.UCSC[Islands.UCSC$Relation_to_Island==region,]
  }
  
  #return X limited to genomic region specified
  return(X[rownames(Islands.UCSC),])
}


###############################################################################################
# average methylation across all probes by a dichotomous variable
###############################################################################################
global.island <- function(beta,pd,split,by='person',full.info=F,non.par=T){
  #beta=beta matrix (row cpgs, col subject)
  #pd=phenotype data
  #binary variable in pd to compare groups
  #if by person: averages for each person, then compares broken up by split
  #if by cpg: average for each CpG, then compares across people broken by split
  
  if(!identical(rownames(pd),colnames(beta))){
    stop('Column names in beta matrix and row names in phenotype data do not match.')
  }
  
  if(by=='person'){
    #average per person
    mean.sub_all <- colMeans(beta)
    mean.sub_sea <- colMeans(genomic.region(beta,'OpenSea'))
    mean.sub_shelf <- colMeans(genomic.region(beta,'Shelf'))
    mean.sub_shore <- colMeans(genomic.region(beta,'Shore'))
    mean.sub_island <- colMeans(genomic.region(beta,'Island'))
    
    #compare
    diff <- data.frame(matrix(0,nrow=5,ncol=6))
    colnames(diff) <- c("feature", "meandiff","upper","lower","p","N")
    diff$feature <- c('All','Open Sea','Shelf','Shore','Island')
    
    pd[,split] <- factor(pd[,split])
    pd[,split] <- relevel(pd[,split],"1")
    
    all <- t.test(mean.sub_all~pd[,split])
    diff[1,2:6] <- c(all$estimate[1]*100-all$estimate[2]*100, all$conf.int[2]*100, all$conf.int[1]*100, all$p.value, length(mean.sub_all))
    
    sea <- t.test(mean.sub_sea~pd[,split])
    diff[2,2:6] <- c(sea$estimate[1]*100-sea$estimate[2]*100, sea$conf.int[2]*100, sea$conf.int[1]*100, sea$p.value, length(mean.sub_sea))
    
    shelf <- t.test(mean.sub_shelf~pd[,split])
    diff[3,2:6] <- c(shelf$estimate[1]*100-shelf$estimate[2]*100, shelf$conf.int[2]*100, shelf$conf.int[1]*100, shelf$p.value, length(mean.sub_shelf))
    
    shore <- t.test(mean.sub_shore~pd[,split])
    diff[4,2:6] <- c(shore$estimate[1]*100-shore$estimate[2]*100, shore$conf.int[2]*100, shore$conf.int[1]*100, shore$p.value, length(mean.sub_shore))
    
    island <- t.test(mean.sub_island~pd[,split])
    diff[5,2:6] <- c(island$estimate[1]*100-island$estimate[2]*100, island$conf.int[2]*100, island$conf.int[1]*100, island$p.value, length(mean.sub_island))
    

    message('Differences are first group - second group:')
    message(paste(names(all$estimate),collapse=' , '))
    
    if(!full.info){
      return(diff)
    }else{
      means <- list(split1=mean.sub_all[pd[,split]==1],split0=mean.sub_all[pd[,split]==0])
      means.island <- list(split1=mean.sub_island[pd[,split]==1],split0=mean.sub_island[pd[,split]==0])
      means.sea <- list(split1=mean.sub_sea[pd[,split]==1],split0=mean.sub_sea[pd[,split]==0])
      return(list(tbl=diff,means=means,means.island=means.island,means.sea=means.sea))
    }
    
    
  }else if(by=='cpg'){
    #split matrix
    beta0 <- beta[,rownames(pd[pd[,split]==0,])]
    beta1 <- beta[,rownames(pd[pd[,split]==1,])]
    
    #average per cpg
    mean.sub_all <- list(rowMeans(beta1),rowMeans(beta0))
    var.sub_all <- list(apply(beta1,1,var),apply(beta0,1,var))
    sea1 <- genomic.region(beta1,'OpenSea')
    sea0 <- genomic.region(beta0,'OpenSea')
    mean.sub_sea <- list(rowMeans(sea1),rowMeans(sea0))
    shelf1 <- genomic.region(beta1, 'Shelf')
    shelf0 <- genomic.region(beta0, 'Shelf')
    mean.sub_shelf <- list(rowMeans(shelf1),rowMeans(shelf0))
    island1 <- genomic.region(beta1,'Island')
    island0 <- genomic.region(beta0,'Island')
    mean.sub_island <- list(rowMeans(island1),rowMeans(island0))
    shore1 <- genomic.region(beta1,'Shore')
    shore0 <- genomic.region(beta0,'Shore')
    mean.sub_shore <- list(rowMeans(shore1),rowMeans(shore0))
    
    #compare
    if(!non.par){
      diff <- data.frame(matrix(0,nrow=5,ncol=6))
      colnames(diff) <- c("feature", "meandiff","upper","lower","p","N")
      diff$feature <- c('All','Open Sea','Shelf','Shore','Island')
      
      all <- t.test(mean.sub_all[[1]],mean.sub_all[[2]],paired=T)
      diff[1,2:6] <- c(all$estimate*100, all$conf.int[2]*100, all$conf.int[1]*100, all$p.value, length(mean.sub_all[[1]]))
      
      sea <- t.test(mean.sub_sea[[1]],mean.sub_sea[[2]],paired=T)
      diff[2,2:6] <- c(sea$estimate*100, sea$conf.int[2]*100, sea$conf.int[1]*100, sea$p.value, length(mean.sub_sea[[1]]))
      
      shelf <- t.test(mean.sub_shelf[[1]],mean.sub_shelf[[2]],paired=T)
      diff[3,2:6] <- c(shelf$estimate*100, shelf$conf.int[2]*100, shelf$conf.int[1]*100, shelf$p.value, length(mean.sub_shelf[[1]]))
      
      shore <- t.test(mean.sub_shore[[1]],mean.sub_shore[[2]],paired=T)
      diff[4,2:6] <- c(shore$estimate*100, shore$conf.int[2]*100, shore$conf.int[1]*100, shore$p.value, length(mean.sub_shore[[1]]))
      
      island <- t.test(mean.sub_island[[1]],mean.sub_island[[2]],paired=T)
      diff[5,2:6] <- c(island$estimate*100, island$conf.int[2]*100, island$conf.int[1]*100, island$p.value, length(mean.sub_island[[1]]))
      
      message('Differences are first group - second group:')
      message(paste(names(all$estimate),collapse=' , '))
      message(paste('N for first group: ',ncol(beta1),sep=''))
      message(paste('N for second group: ',ncol(beta0),sep=''))
      
      if(!full.info){
        return(diff)
      }else{
        names(mean.sub_all) <- c('split1','split0')
        names(mean.sub_island) <- c('split1','split0')
        names(mean.sub_sea) <- c('split1','split0')
        return(list(tbl=diff,means=mean.sub_all,means.island=mean.sub_island,means.sea=mean.sub_sea))
      }
    }else{
      diff <- data.frame(matrix(0,nrow=5,ncol=3))
      colnames(diff) <- c("feature","p","N")
      diff$feature <- c('All','Open Sea','Shelf','Shore','Island')
      
      all <- wilcox.test(mean.sub_all[[1]],mean.sub_all[[2]])
      diff[1,2:3] <- c(all$p.value, length(mean.sub_all[[1]]))
      
      sea <- wilcox.test(mean.sub_sea[[1]],mean.sub_sea[[2]])
      diff[2,2:3] <- c(sea$p.value, length(mean.sub_sea[[1]]))
      
      shelf <- wilcox.test(mean.sub_shelf[[1]],mean.sub_shelf[[2]])
      diff[3,2:3] <- c(shelf$p.value, length(mean.sub_shelf[[1]]))
      
      shore <- wilcox.test(mean.sub_shore[[1]],mean.sub_shore[[2]])
      diff[4,2:3] <- c(shore$p.value, length(mean.sub_shore[[1]]))
      
      island <- wilcox.test(mean.sub_island[[1]],mean.sub_island[[2]])
      diff[5,2:3] <- c(island$p.value, length(mean.sub_island[[1]]))
      
      message('Non-parametric test: no direction of effects')
      
      return(diff)
    }
  }else if(by=='cpg.secret'){
    #split matrix
    beta0 <- beta[,rownames(pd[pd[,split]==0,])]
    beta1 <- beta[,rownames(pd[pd[,split]==1,])]
    
    #average per cpg
    mean.sub_all <- list(rowMeans(beta0),rowMeans(beta1))
    mean.sub_sea <- list(rowMeans(genomic.region(beta0,'OpenSea')),rowMeans(genomic.region(beta1,'OpenSea')))
    mean.sub_nshelf <- list(rowMeans(genomic.region(beta0,'N_Shelf')),rowMeans(genomic.region(beta1,'N_Shelf')))
    mean.sub_nshore <- list(rowMeans(genomic.region(beta0,'N_Shore')),rowMeans(genomic.region(beta1,'N_Shore')))
    mean.sub_island <- list(rowMeans(genomic.region(beta0,'Island')),rowMeans(genomic.region(beta1,'Island')))
    mean.sub_sshore <- list(rowMeans(genomic.region(beta0,'S_Shore')),rowMeans(genomic.region(beta1,'S_Shore')))
    mean.sub_sshelf <- list(rowMeans(genomic.region(beta0,'S_Shelf')),rowMeans(genomic.region(beta1,'S_Shelf')))
    
    #compare
    diff <- data.frame(matrix(0,nrow=7,ncol=6))
    colnames(diff) <- c("feature", "meandiff","upper","lower","p","N")
    diff$feature <- c('All','Open Sea','North Shelf','North Shore','Island','South Shore','South Shelf')
    
    all <- t.test(mean.sub_all[[1]],mean.sub_all[[2]])
    diff[1,2:6] <- c(all$estimate[1]*100-all$estimate[2]*100, all$conf.int[2]*100, all$conf.int[1]*100, all$p.value, length(mean.sub_all[[1]]))
    
    sea <- t.test(mean.sub_sea[[1]],mean.sub_sea[[2]])
    diff[2,2:6] <- c(sea$estimate[1]*100-sea$estimate[2]*100, sea$conf.int[2]*100, sea$conf.int[1]*100, sea$p.value, length(mean.sub_sea[[1]]))
    
    nshelf <- t.test(mean.sub_nshelf[[1]],mean.sub_nshelf[[2]])
    diff[3,2:6] <- c(nshelf$estimate[1]*100-nshelf$estimate[2]*100, nshelf$conf.int[2]*100, nshelf$conf.int[1]*100, nshelf$p.value, length(mean.sub_nshelf[[1]]))
    
    nshore <- t.test(mean.sub_nshore[[1]],mean.sub_nshore[[2]])
    diff[4,2:6] <- c(nshore$estimate[1]*100-nshore$estimate[2]*100, nshore$conf.int[2]*100, nshore$conf.int[1]*100, nshore$p.value, length(mean.sub_nshore[[1]]))
    
    island <- t.test(mean.sub_island[[1]],mean.sub_island[[2]])
    diff[5,2:6] <- c(island$estimate[1]*100-island$estimate[2]*100, island$conf.int[2]*100, island$conf.int[1]*100, island$p.value, length(mean.sub_island[[1]]))
    
    sshore <- t.test(mean.sub_sshore[[1]],mean.sub_sshore[[2]])
    diff[6,2:6] <- c(sshore$estimate[1]*100-sshore$estimate[2]*100, sshore$conf.int[2]*100, sshore$conf.int[1]*100, sshore$p.value, length(mean.sub_sshore[[1]]))
    
    sshelf <- t.test(mean.sub_sshelf[[1]],mean.sub_sshelf[[2]])
    diff[7,2:6] <- c(sshelf$estimate[1]*100-sshelf$estimate[2]*100, sshelf$conf.int[2]*100, sshelf$conf.int[1]*100, sshelf$p.value, length(mean.sub_sshelf[[1]]))
    
    message('Differences are first group - second group:')
    message(paste(names(all$estimate),collapse=' , '))
    message(paste('N for first group: ',ncol(beta0),sep=''))
    message(paste('N for second group: ',ncol(beta1),sep=''))
    
    return(diff)
  }else {
    stop('Please choose person or cpg for by parameter')
  }
  
}

###############################################################################################
# remove probes in sex chromosomes
###############################################################################################
autosomes <- function(X, anno='450k'){
  #X - something with cpg as rownames
  #anno = 450k or epic
  
  if(anno=='450k' | anno=='450K'){
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }else if(anno=='epic' | anno=='EPIC'){
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  }else{
    stop('Dunno that array m8, 450k or EPIC only.')
  }
  
  data(Locations)
  Locations <- Locations[Locations$chr=='chrX' | Locations$chr=='chrY',]
  
  return(X[!rownames(X)%in%rownames(Locations),])
}

###############################################################################################
# plot density and histograms
###############################################################################################
plot.denhist <- function(global, file='NO_NAME_GIVEN.pdf', means='means'){
  pdf(file)
  split1 <- global[[means]][['split1']]
  split1 <- cbind(rep('ASD Case',length(split1)),split1)
  split0 <- global[[means]][['split0']]
  split0 <- cbind(rep('Non-ASD',length(split0)),split0)
  means <- rbind(split1,split0)
  rownames(means) <- NULL
  means <- data.frame(means)
  names(means) <- c('case','mean')
  means$mean <- as.numeric(as.character(means$mean))
  print(ggplot(means, aes(mean, fill=case)) + geom_density(alpha=0.3))
  print(ggplot(means, aes(mean, fill=case)) + geom_histogram(alpha=0.3, aes(), position='identity'))
  dev.off()
  'Check those plots yo'
}


#beta = beta matrix
#pd = phenotype file
#split = variable of interest for bivariate analysis
bivariate.global <- function(beta, pd, split, med=F){
  
  beta1 <- beta[,colnames(beta) %in% rownames(pd[pd[,split]==1,])]
  beta0 <- beta[,colnames(beta) %in% rownames(pd[pd[,split]==0,])]
  col1 <- ncol(beta1)
  col0 <- ncol(beta0)
  
  #cpg annotation
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  an450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  island_enh <- an450k[,c('Relation_to_Island', 'Enhancer')]
  beta1 <- merge(beta1,island_enh,by='row.names',sort=F,all.x=T,all.y=F)
  beta0 <- merge(beta0,island_enh,by='row.names',sort=F,all.x=T,all.y=F)
  if(!med){
    beta1$Rowmeans <- rowMeans(beta1[,2:(col1+1)])
    beta0$Rowmeans <- rowMeans(beta0[,2:(col0+1)])
  }else{
    library(matrixStats)
    beta1$Rowmeans <- rowMedians(as.matrix(beta1[,2:(col1+1)]))
    beta0$Rowmeans <- rowMedians(as.matrix(beta0[,2:(col0+1)]))
  }
  
  
  return(list(beta1,beta0))
}

#global methylation by relation of islands
island.relation <- function(results) {
  res <- data.frame(matrix(NA, nrow=8, ncol=6))
  colnames(res) <- c("feature", "meandiff","upper","lower","p","N")
  res$feature <- c("All", "Island","North Shore","South Shore","North Shelf","South Shelf","Open Sea","Enhancer")
  test <- t.test(results$logFC)
  res[1,2]<-test$est*100
  res[1,3]<-test$conf.int[1]*100
  res[1,4]<-test$conf.int[2]*100
  res[1,5]<-test$p.value[1]
  test<-t.test(results$logFC[results$island=="Island"])
  res[2,2]<-test$est*100
  res[2,3]<-test$conf.int[1]*100
  res[2,4]<-test$conf.int[2]*100
  res[2,5]<-test$p.value[1]
  test<-t.test(results$logFC[results$island=="N_Shore"])
  res[3,2]<-test$est*100
  res[3,3]<-test$conf.int[1]*100
  res[3,4]<-test$conf.int[2]*100
  res[3,5]<-test$p.value[1]
  test<-t.test(results$logFC[results$island=="S_Shore"])
  res[4,2]<-test$est*100
  res[4,3]<-test$conf.int[1]*100
  res[4,4]<-test$conf.int[2]*100
  res[4,5]<-test$p.value[1]
  test<-t.test(results$logFC[results$island=="N_Shelf"])
  res[5,2]<-test$est*100
  res[5,3]<-test$conf.int[1]*100
  res[5,4]<-test$conf.int[2]*100
  res[5,5]<-test$p.value[1]
  test<-t.test(results$logFC[results$island=="S_Shelf"])
  res[6,2]<-test$est*100
  res[6,3]<-test$conf.int[1]*100
  res[6,4]<-test$conf.int[2]*100
  res[6,5]<-test$p.value[1]
  test<-t.test(results$logFC[results$island=="OpenSea"])
  res[7,2]<-test$est*100
  res[7,3]<-test$conf.int[1]*100
  res[7,4]<-test$conf.int[2]*100
  res[7,5]<-test$p.value[1]
  test<-t.test(results$logFC[results$enhancer==TRUE])
  res[8,2]<-test$est*100
  res[8,3]<-test$conf.int[1]*100
  res[8,4]<-test$conf.int[2]*100
  res[8,5]<-test$p.value[1]
  res$N<-c(length(results$logFC), length(results$logFC[results$island=="Island"]), length(results$logFC[results$island=="N_Shore"]), length(results$logFC[results$island=="S_Shore"]), length(results$logFC[results$island=="N_Shelf"]), length(results$logFC[results$island=="S_Shelf"]), length(results$logFC[results$island=="OpenSea"]), length(results$logFC[results$enhancer==TRUE]))
  
  res$p <- ifelse(res$p==0, 2.2e-16, res$p)
  
  return(res[c(1,7,6,4,2,3,5,8),])
}

#global methylation from beta matrix split by categorical variable
island.relation.matrix <- function(beta1, beta2) {
  res <- data.frame(matrix(NA, nrow=8, ncol=6))
  colnames(res) <- c("feature", "meandiff","upper","lower","p","N")
  res$feature <- c("All", "Island","North Shore","South Shore","North Shelf","South Shelf","Open Sea","Enhancer")
  test <- t.test(beta1$Rowmeans, beta2$Rowmeans, paired=T)
  res[1,2]<-test$est*100
  res[1,3]<-test$conf.int[1]*100
  res[1,4]<-test$conf.int[2]*100
  res[1,5]<-test$p.value[1]
  test<-t.test(beta1$Rowmeans[beta1$Relation_to_Island=="Island"], beta2$Rowmeans[beta2$Relation_to_Island=="Island"], paired=T)
  res[2,2]<-test$est*100
  res[2,3]<-test$conf.int[1]*100
  res[2,4]<-test$conf.int[2]*100
  res[2,5]<-test$p.value[1]
  test<-t.test(beta1$Rowmeans[beta1$Relation_to_Island=="N_Shore"], beta2$Rowmeans[beta2$Relation_to_Island=="N_Shore"], paired=T)
  res[3,2]<-test$est*100
  res[3,3]<-test$conf.int[1]*100
  res[3,4]<-test$conf.int[2]*100
  res[3,5]<-test$p.value[1]
  test<-t.test(beta1$Rowmeans[beta1$Relation_to_Island=="S_Shore"], beta2$Rowmeans[beta2$Relation_to_Island=="S_Shore"], paired=T)
  res[4,2]<-test$est*100
  res[4,3]<-test$conf.int[1]*100
  res[4,4]<-test$conf.int[2]*100
  res[4,5]<-test$p.value[1]
  test<-t.test(beta1$Rowmeans[beta1$Relation_to_Island=="N_Shelf"], beta2$Rowmeans[beta2$Relation_to_Island=="N_Shelf"], paired=T)
  res[5,2]<-test$est*100
  res[5,3]<-test$conf.int[1]*100
  res[5,4]<-test$conf.int[2]*100
  res[5,5]<-test$p.value[1]
  test<-t.test(beta1$Rowmeans[beta1$Relation_to_Island=="S_Shelf"], beta2$Rowmeans[beta2$Relation_to_Island=="S_Shelf"], paired=T)
  res[6,2]<-test$est*100
  res[6,3]<-test$conf.int[1]*100
  res[6,4]<-test$conf.int[2]*100
  res[6,5]<-test$p.value[1]
  test<-t.test(beta1$Rowmeans[beta1$Relation_to_Island=="OpenSea"], beta2$Rowmeans[beta2$Relation_to_Island=="OpenSea"], paired=T)
  res[7,2]<-test$est*100
  res[7,3]<-test$conf.int[1]*100
  res[7,4]<-test$conf.int[2]*100
  res[7,5]<-test$p.value[1]
  test<-t.test(beta1$Rowmeans[beta1$Enhancer==TRUE], beta2$Rowmeans[beta2$Enhancer==TRUE], paired=T)
  res[8,2]<-test$est*100
  res[8,3]<-test$conf.int[1]*100
  res[8,4]<-test$conf.int[2]*100
  res[8,5]<-test$p.value[1]
  res$N<-c(length(beta1$Rowmeans), length(beta1$Rowmeans[beta1$Relation_to_Island=="Island"]), length(beta1$Rowmeans[beta1$Relation_to_Island=="N_Shore"]), length(beta1$Rowmeans[beta1$Relation_to_Island=="S_Shore"]), length(beta1$Rowmeans[beta1$Relation_to_Island=="N_Shelf"]), length(beta1$Rowmeans[beta1$Relation_to_Island=="S_Shelf"]), length(beta1$Rowmeans[beta1$Relation_to_Island=="OpenSea"]), length(beta1$Rowmeans[beta1$Enhancer==TRUE]))
  
  res$p <- ifelse(res$p==0, 2.2e-16, res$p)
  
  return(res[c(1,7,6,4,2,3,5,8),])
}

###################################################
#gene enrichment with arbitrary set with missMethyl
###################################################
gset.enrich <- function(res,gset,sparse.p=T,missMethyl=T){
  if(missMethyl){
  library(missMethyl)
  
  #set p value cutoffs for results
  if(sparse.p){
    cutoffs <- c(0.0001, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04 , 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.99)
  }
  else{
    cutoffs <- c(seq(-4,-2,0.1),seq(-1.95,-0,0.05))
    cutoffs <- 10^cutoffs
    cutoffs[length(cutoffs)] <- 0.99
  }
  
  #prepare blank results to be filled in
  enrich <- as.list(rep(0,length(cutoffs)))
  names(enrich) <- paste('EWAS P ',cutoffs,sep='')
  p.values <- rep(0,length(cutoffs))
  
  #run enrichment test for each cutoff
  for(i in 1:length(cutoffs)){
    print(paste0('Now doing cutoff P ',cutoffs[i]))
    sig <- rownames(res[res$P.Value<cutoffs[i],])
    all <- rownames(res)
    test <- gsameth(sig,all,collection=gset,array.type='450K')
    print(paste0('-------------------------'))
    enrich[[i]] <- test
    p.values[i] <- -log(test[,'P.DE'],10)
  }
  
  return(list(tables=enrich,cutoffs=cutoffs,logp=p.values))
  }else
  {
    
  ### same as in the sfari.enrich method, fisher's test
  annot<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/NCSv/NoControls/Resid/annot-loc-symbol-entrez09-24-13")
  dat.annot<-merge(res, annot, by.x="row.names",by.y="X", sort=F)
  
  #First take lowest P value per gene
  dat.na<-dat.annot[!(is.na(dat.annot$entrezids)),]
  dat.o<-dat.na[order(dat.na$entrezids, dat.na$P.Value),]
  dat.o$index<-ave(rep(NA, nrow(dat.o)), dat.o$entrezids, FUN=seq_along) 
  dat.low.p<-dat.o[dat.o$index==1,]
  
  dat.low.p$gset<-ifelse(dat.low.p$symbols %in% gset, 1, 0)
  
  library(MASS)  
  
  if(sparse.p){
    cutoffs <- c(0.0001, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04 , 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.99)
  }
  else{
    cutoffs <- c(seq(-4,-2,0.1),seq(-1.95,-0,0.05))
    cutoffs <- 10^cutoffs
    cutoffs[length(cutoffs)] <- 0.99
  }
  enrich <- as.list(rep(0,length(cutoffs)))
  fish <- as.list(rep(0,length(cutoffs)))
  names(enrich) <- paste('EWAS P ',cutoffs,sep='')
  p.values <- rep(0,length(cutoffs))
  
  for(i in 1:length(cutoffs)){
    pcut <- ifelse(dat.low.p$P.Value<cutoffs[i],1,0)
    tbl <- table(pcut,gset = dat.low.p$gset)
    test <- fisher.test(tbl)
    enrich[[i]] <- tbl
    fish[[i]] <- test$estimate
    p.values[i] <- -log(test$p.value,10)
  }
  
  return(list(tables = enrich, cutoffs = cutoffs, logp = p.values))
  }
}

#dat = single site results from limma regression
#returns list with 2x2 tables, p-value cutoffs and enrichment significance
#gene.overlap.at.peak: if true, returns list of genes in overlap with sfari at peak significance cutoff
#fisher: if true, uses fisher exact test
#sparse.p: if true, uses a smaller list of p-value cutoffs
#gene.all: if true, includes list of all genes with CpG at cut
sfari.enrich <- function(dat, gene.overlap.at.peak = T, fisher=T, sparse.p=T, cut=NULL,gene.at.cut=T){
  #sfari gene file
  #sfari<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/SFARI/SFARI-gene-summary-20150123.csv")
  sfari<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/SFARI-Gene_genes_export28-08-2017.csv")
  sfari.na<-sfari[!(is.na(sfari$gene.symbol)),]
  
  #Illumina product annotation (includes unannotated probes)
  annot<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/NCSv/NoControls/Resid/annot-loc-symbol-entrez09-24-13")
  dat.annot<-merge(dat, annot, by.x="row.names",by.y="X", sort=F)
  
  #   #annotation nearest gene (all have an annotation)
  #   neargene<-read.csv("/dcl01/NDEpi/data/Projects/InProgress/kbakulsk/kbakulsk/Illumina450k-annotateneargene-20150311.csv")
  #   #neargene<-read.csv("/Users/kbakulsk/Google Drive/Env_DNAm_Enrichment/Illumina450k-annotateneargene-20150311.csv")
  #   dat.near<-merge(dat, neargene, by="X", sort=F)
  
  #First take lowest P value per gene
  dat.na<-dat.annot[!(is.na(dat.annot$entrezids)),]
  dat.o<-dat.na[order(dat.na$entrezids, dat.na$P.Value),]
  dat.o$index<-ave(rep(NA, nrow(dat.o)), dat.o$entrezids, FUN=seq_along) 
  dat.low.p<-dat.o[dat.o$index==1,]
  
  dat.low.p$sfari<-ifelse(dat.low.p$symbols %in% sfari$gene.symbol, 1, 0)

  
  library(MASS)  
  
  if(sparse.p){
    cutoffs <- c(0.0001, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04 , 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.99)
  }
  else{
    cutoffs <- c(seq(-4,-2,0.1),seq(-1.95,-0,0.05))
    cutoffs <- 10^cutoffs
    cutoffs[length(cutoffs)] <- 0.99
  }
  enrich <- as.list(rep(0,length(cutoffs)))
  fish <- as.list(rep(0,length(cutoffs)))
  names(enrich) <- paste('EWAS P ',cutoffs,sep='')
  p.values <- rep(0,length(cutoffs))
  
  for(i in 1:length(cutoffs)){
    pcut <- ifelse(dat.low.p$P.Value<cutoffs[i],1,0)
    tbl <- table(pcut,sfari = dat.low.p$sfari)
    if(fisher){test <- fisher.test(tbl)
    } else {test <- chisq.test(tbl)}
    enrich[[i]] <- tbl
    fish[[i]] <- test$estimate
    p.values[i] <- -log(test$p.value,10)
  }
  

  if(gene.overlap.at.peak){
    gene.list <- sfari.whichgene(dat.low.p,p.values,cutoffs,cut)
    expected <- chisq.test(enrich[[paste('EWAS P ',gene.list$maxP,sep='')]])$expected
    if(!gene.at.cut){
      return(list(tables = enrich, OR = fish, cutoffs = cutoffs, p.values = p.values,
                genes = gene.list$genes, peak.p.cutoff = gene.list$maxP, expected = expected))
    }else{
      dat.low.p$pcut <- ifelse(dat.low.p$P.Value<cut,1,0)
      genes.at.cut <- as.character(dat.low.p[dat.low.p$pcut==1,'symbols'])
      return(list(tables = enrich, OR = fish, cutoffs = cutoffs, p.values = p.values,
                genes = gene.list$genes, peak.p.cutoff = gene.list$maxP, expected = expected,genes.at.cut=genes.at.cut))
    }
  }
  
  return(list(tables = enrich, OR = fish, cutoffs = cutoffs, p.values = p.values))
}

#gives which genes are overlapped with sfari list at most significant p cut-off
#unless a specific cut point is given to override
#used by the sfari.enrich function
sfari.whichgene <- function(dat.low.p, p.values, cutoffs, cut){
  if(is.null(cut)){
    maxIndex <- which(p.values==max(p.values))
    maxP <- cutoffs[maxIndex]
  } else {
    maxP <- cut
  }
  dat.low.p$pcut <- ifelse(dat.low.p$P.Value<maxP,1,0)
  genes <- dat.low.p[dat.low.p$pcut==1 & dat.low.p$sfari==1,]
  genes <- genes[,c('entrezids','symbols')]
  return(list(genes=genes,maxP=maxP))
}