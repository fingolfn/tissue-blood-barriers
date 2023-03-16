
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########      ########      ########     ########             #########           ##       ############      ##
##     ##     ##     ##     ##     ##           ##     ##            ##      ##         ####           ##          ####
##     ##     ##     ##    ##      ##           ##     ##            ##       ##       ##  ##          ##         ##  ##
##     ########      #######       ########     ########             ##        ##     ##    ##         ##        ##    ##
##     ##            ##    ##      ##           ##                   ##       ##     ##########        ##       ##########
##     ##            ##     ##     ##           ##                   ##      ##     ##        ##       ##      ##        ##
##     ##            ##      ##    ########     ##                   #########     ##          ##      ##     ##          ##
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#GOAL: Process "References Dataset" into:
#         -> BR-matrix (for App:  https://douglasslab.shinyapps.io/barrier_ratio/)
#         -> KO-matrix (for App: https://douglasslab.shinyapps.io/ko_ratio/)
#         -> BBB-matrix (for figure 3 analysis) 

##############################################################################################
#1   reformat raw data (so can easily populate BR and KO matrices)
##############################################################################################

#LOAD RAW DATA:
tmp_raw<-read.csv("raw_data/REFERENCES-meta-review-01_text.csv",as.is=T,header=T)


#RELABEL column names to make easier to work with:
tissue2label<-c("brain","muscle","colon","caecum","intestines","stomach","liver","kidney","lungs",
                "spleen","heart","ovary","uterus","breast","testes","gall","fat.neck","fat.organ","epododumis","thymus","lymph","penis","brain.Cu")

names(tissue2label)<-c("Brain","Muscle","Colon","Caecum","Small Gut","Stomach","Liver","Kidney","Lung",
                       "Spleen","Heart","Ovary","Uterus","Breast","Teste","Gall Bladder","Fat (neck)","Fat (organ)","Epododumis","Thymus","Lymph Node","Penis","Brain - unbound")

tmp_raw$TISSUE<-tissue2label[tmp_raw$TISSUE]


#RELABEL:  1a/1b ko -> 1a ko (to keep annotation simple given 1b not appreciably expressed in tissue-barriers (mostly stem cells))
tmp_raw$GENE<-gsub("abcb1a/1b -/-","abcb1a -/-",tmp_raw$GENE)


##############################################################################################
#2 MAKE BR MATRIX:
##############################################################################################

#extract data for top 4 genes
genes_sorted<-sort(table(tmp_raw$GENE),decreasing=T)
tmp_raw_topGenes<-tmp_raw[which(tmp_raw$GENE %in% names(genes_sorted[1:4])),]

#Make unique experiment-names column (will be rownames for BR-matrix)
tmp_raw_topGenes$BR_names<-paste(tmp_raw_topGenes$DRUG,
                                 tmp_raw_topGenes$Method,
                                 tmp_raw_topGenes$dose..mg.kg.,
                                 tmp_raw_topGenes$Admin.Route,
                                 tmp_raw_topGenes$time,
                                 tmp_raw_topGenes$GENE,sep="_")



#make BR-number matrix for each experiment
unique_expr<-unique(tmp_raw_topGenes$BR_names)

BR_matrix_numbers<-matrix(NA,ncol=length(tissue2label),nrow=length(unique_expr))
rownames(BR_matrix_numbers)<-unique_expr
colnames(BR_matrix_numbers)<-tissue2label


for (row in rownames(BR_matrix_numbers)){
  for (col in colnames(BR_matrix_numbers)){
    br_idx<-which(tmp_raw_topGenes$BR_names==row & tmp_raw_topGenes$TISSUE==col)
    
    if (length(br_idx)==1){
      BR_matrix_numbers[row,col]<-tmp_raw_topGenes$Barrier.ratio[br_idx]
    }
    
    #added this as ponatinib had replicates w/ same stats
    if (length(br_idx)>1){
      BR_matrix_numbers[row,col]<-tmp_raw_topGenes$Barrier.ratio[br_idx[1]]
      print(row)
    }
    
    
  }
}

#make Experimental-details matrix
BR_matrix_text<-cbind(sapply(strsplit(rownames(BR_matrix_numbers),split="_"), function(x) x[1]),
                      sapply(strsplit(rownames(BR_matrix_numbers),split="_"), function(x) x[3]),
                      sapply(strsplit(rownames(BR_matrix_numbers),split="_"), function(x) x[4]),
                      sapply(strsplit(rownames(BR_matrix_numbers),split="_"), function(x) x[6]),
                      sapply(strsplit(rownames(BR_matrix_numbers),split="_"), function(x) x[5]),
                      sapply(strsplit(rownames(BR_matrix_numbers),split="_"), function(x) x[2]))
rownames(BR_matrix_text)<-rownames(BR_matrix_numbers)
colnames(BR_matrix_text)<-c("drug","dose","admin","gene","time","method")

#combine BR-number and Experimental-details
BR_matrix<-cbind(BR_matrix_text[rownames(BR_matrix_numbers),],BR_matrix_numbers)


#SAVE DATA:  R and csv (for excel for supp info)
save(BR_matrix,file="BR_matrix.RData")

#write.csv(BR_matrix,file="BR_matrix.csv")






##############################################################################################
#3   MAKE KO MATRIX:
##############################################################################################

#LOAD BR-matrix (made in step #2)
load("BR_matrix.RData")

#SPLIT BR-matrix into WT-matrix and KO-matrix to divide into KO-matrix
WT_matrix<-BR_matrix[which(BR_matrix[,"gene"]=="WT"),]
KO_matrix<-BR_matrix



#row_idx=1
for (row_idx in 1:nrow(KO_matrix)){
  ko_drug=KO_matrix[row_idx,"drug"]
  ko_dose=KO_matrix[row_idx,"dose"]
  ko_admin=KO_matrix[row_idx,"admin"]
  ko_time=KO_matrix[row_idx,"time"]
  ko_method=KO_matrix[row_idx,"method"]
  
  wt_idx<-which(WT_matrix[,"drug"]==ko_drug &  WT_matrix[,"time"]==ko_time &  WT_matrix[,"method"]==ko_method &  WT_matrix[,"dose"]==ko_dose &  WT_matrix[,"admin"]==ko_admin)

  
  if(length(wt_idx)>0){
    #if multiple WT measurements for drug & time
    if(length(wt_idx)>1){
      wt_idx<-wt_idx[1]
      print(ko_drug)
    }
    
    KO_matrix[row_idx,7:29]<-as.numeric(BR_matrix[row_idx,7:29])/as.numeric(WT_matrix[wt_idx,7:29])
  }else{
    KO_matrix[row_idx,7:29]<-rep(NA,length(7:29))
  }
  
}

#REMOVE ROWS THAT DIDN'T HAVE WT REF:
KO_matrix<-KO_matrix[-which(is.na(KO_matrix[,"brain"])==T),]


#SAVE DATA:
save(KO_matrix,file="KO_matrix.RData")

#write.csv(KO_matrix,file="KO_matrix.csv")





##############################################################################################
#3   MAKE BBB MATRIX:
##############################################################################################


load("BR_matrix.RData")



#no NA's in BBB matrix:
bbb_na_idx<-which(is.na(BR_matrix[,"brain"])==T)


drug_time_inorder<-paste(BR_matrix[,"drug"],paste0(BR_matrix[,"time"],"hr"),sep="_")


BBB_matrix<-matrix(NA,nrow=length(unique(drug_time_inorder)),ncol=length(unique(BR_matrix[,"gene"])))
rownames(BBB_matrix)<-unique(drug_time_inorder)
colnames(BBB_matrix)<-unique(BR_matrix[,"gene"])[c(4,1:3)]

#row="Doxorubicin_24hr"
for (row in rownames(BBB_matrix)){
  #col="WT"
  for (col in colnames(BBB_matrix)){
    data_idx<-which(drug_time_inorder==row & BR_matrix[,"gene"]==col)
    
    if (length(data_idx)==1){
      BBB_matrix[row,col]<-as.numeric(BR_matrix[data_idx,"brain"])
    }
    
    if (length(data_idx)>1){
      data_idx=data_idx[1]
      BBB_matrix[row,col]<-as.numeric(BR_matrix[data_idx,"brain"])
    }
    
  }
}



save(BBB_matrix,file="BBB_matrix.RData")

#write.csv(BBB_matrix,file="BBB_matrix.csv")









####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########      ########                         #########           ##       ############      ##                            
##     ##     ##     ##     ##                        ##      ##         ####           ##          ####
##     ##     ##     ##    ##                         ##       ##       ##  ##          ##         ##  ##
##     ########      #######          #########       ##        ##     ##    ##         ##        ##    ##
##     ##     ##     ##    ##                         ##       ##     ##########        ##       ##########
##     ##     ##     ##     ##                        ##      ##     ##        ##       ##      ##        ##
##     ########      ##      ##                       ########      ##          ##      ##     ##          ##
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ##      ###     ##      ##      ##       ##      ##     #######     ##    ########                                 
##      ####     ####    ##     ####     ##        ##    ##     ##           ##   ##                         
##     ##  ##    ## ##   ##    ##  ##    ##         ##  ##      ##           ##   ##                           
##    ##    ##   ##  ##  ##   ##    ##   ##          ####        #######     ##    #######                                         
##    ########   ##   ## ##   ########   ##           ##               ##    ##          ##                            
##    ##    ##   ##    ####   ##    ##   ##           ##               ##    ##          ##                   
##    ##    ##   ##     ###   ##    ##   #######      ##         #######     ##    #######                                              
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#GOAL:  analysis on raw data

load("BR_matrix.RData")
load("KO_matrix.RData")
load("BBB_matrix.RData")


#NEEDED LIBRARIES:
library("gplots")
library("RColorBrewer")

#PLOTTING FUNCTION:
hcluster_unsortBW<-function(matrix,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="binary"){
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=NA,Colv=NA, cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("white", "black"))(n = 2),
            density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab = units,symm=F,dendrogram="none")
  
  return(rowMeans(matrix))
}
hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "log10(ratio)",key.xlab = "log10(ratio)")
}

#LINEAR REGRESSION:
linear_fit_plot<-function(x,y,pt_colors="black",xlab=NA,ylab=NA){
  #linear regression:
  lin_reg<-lm(y~x)
  #make fit equation + r-squared:
  pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
  fit_coeff<-round(coef(lin_reg),6)
  r2 <- round(summary(lin_reg)$r.squared, 2)
  rmse <- round(sqrt(mean(resid(lin_reg)^2)), 2)
  eq_r2<-paste("r = ", r2," "," "," ","p-value = ",signif(pvalue_corr,2))
  #plot data:
  plot(x,y, pch = 16, cex = 0.8, col = pt_colors,xlab=xlab,ylab=ylab)
  #add fit-line and fit-equation:
  abline(lin_reg)
  mtext(eq_r2, 3, line=-2,cex=1)
}







##############################################################################################################
#FIGURE 2:   understanding biases in dataset
##############################################################################################################

#############################
#Visualize gaps in data
#############################

#EXTRACT DATA:
BR_matrix_data<-data.matrix(BR_matrix[,7:29])
class(BR_matrix_data)<-"numeric"

#Make Gaps Data (replace NA w/ 0 and non-NA with 1)
BR_matrix_gaps<-BR_matrix_data
BR_matrix_gaps[which(is.na(BR_matrix_data)==T)]<-0
BR_matrix_gaps[which(is.na(BR_matrix_data)==F)]<-1

#sort rows and columns by number of NA's
rows_sort<-names(sort(rowSums(BR_matrix_gaps),decreasing=T))
cols_sort<-names(sort(colSums(BR_matrix_gaps),decreasing=T))

hcluster_unsortBW(t(BR_matrix_gaps[rows_sort,cols_sort]),col_size = 0.2,title="Measurements from 422 Expr")

#Total Data available
sum(BR_matrix_gaps)

#############################
#Figure 2B:   Total Measurements/Tissue
#############################
par(las=2,mar=c(5, 7, 4, 2)) # make label text perpendicular to axis
barplot(sort(colSums(BR_matrix_gaps),decreasing=F), main="Measurements/Tissue", beside=F,cex.axis = 1,cex=1,cex.main=2,col="black",horiz=T,log="x",xlim=c(1,500))



#############################
#Figure 2C total measurements / drug (top 20 drugs)
#############################
measure_per_expr<-sort(rowSums(BR_matrix_gaps),decreasing=T)
names(measure_per_expr)<-sapply(strsplit(names(measure_per_expr),split="_"),function(x) x[1])
expr_per_drug<-sort(tapply(measure_per_expr,names(measure_per_expr),sum))

par(las=2,mar=c(5, 7, 4, 2)) # make label text perpendicular to axis
barplot(expr_per_drug[53:73], main="Total Measurements/Drug", beside=F,cex.axis = 1,cex=1,
        cex.main=2,col="black",horiz=T,log="x",xlim=c(8,100))





#############################
#total measurements / gene
#############################
measureGene_per_expr<-sort(rowSums(BR_matrix_gaps),decreasing=T)
names(measureGene_per_expr)<-sapply(strsplit(names(measureGene_per_expr),split="_"),function(x) x[6])
expr_per_gene<-sort(tapply(measureGene_per_expr,names(measureGene_per_expr),sum))

par(las=2,mar=c(5, 7, 4, 2)) # make label text perpendicular to axis
barplot(expr_per_gene, main="Total Measurements/Gene", beside=F,cex.axis = 1,cex=1,
        cex.main=2,col="black",horiz=T,log="x",xlim=c(100,350))




##############################################################################################################
#FIGURE 3:   Blood brain barrier:   ~40 drugs x 4 mouse models (WT, B1, G2, TKO) 
##############################################################################################################
#need to make BBB matrix for figure


#make gaps matrix
BBB_matrix_gaps<-BBB_matrix
BBB_matrix_gaps[which(is.na(BBB_matrix)==T)]<-0
BBB_matrix_gaps[which(is.na(BBB_matrix)==F)]<-1

#ID rows with no NA's
rows_all_data<-names(which(sort(rowSums(BBB_matrix_gaps),decreasing=T)==4))

#Rank drugs by BBB-ratio
ranked_bbb_drugs<-names(sort(BBB_matrix[rows_all_data,"TKO"]/BBB_matrix[rows_all_data,"WT"],decreasing=T))

#SET MAX (100) & MIN (0.01) for consistent color-scales
BBB_matrix_tmp<-BBB_matrix
BBB_matrix_tmp[which(BBB_matrix<(0.01))]<-0.01
BBB_matrix_tmp[which(BBB_matrix>(100))]<-100

#PLOT:
dev.off()
hcluster(log10(BBB_matrix_tmp[ranked_bbb_drugs[1:40],]),row_size = 0.7)





##############################################################################################################
#FIGURE 5:   B1 Tissue Redistribution:  top drugs tissues
##############################################################################################################
#ID rows for Abcb1a -/- data for top drugs
rel_idx<-which(KO_matrix[,"gene"]=="abcb1a -/-" & KO_matrix[,"drug"] %in% c("Cyclosporin A","Doxorubicin","Dexamethasone","Digoxin","Ivermectin","Morphine","Vinblastine"))

#EXTRACT KO data:      
KO_matrix_data<-data.matrix(KO_matrix[rel_idx,7:29])
class(KO_matrix_data)<-"numeric"


#SET MAX/MIN for color-scale
KO_matrix_data[which(KO_matrix_data<(0.01))]<-0.01
KO_matrix_data[which(KO_matrix_data>(100))]<-100

KO_matrix_data_invert<-(KO_matrix_data)^-1

hcluster(log10(KO_matrix_data_invert[,c("brain","muscle","intestines","stomach","liver","kidney","lungs","spleen","heart","testes","gall","thymus")]))





####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########     ###     ##         ##                                   
##     ##     ##    ####    ##        ####            
##     ##    ##     ## ##   ##       ##  ##       #####  #####    #####
##     #######      ##  ##  ##      ##    ##     ##      ##      ##   ##
##     ##    ##     ##   ## ##     ##########     ####   #####   ##  ###
##     ##     ##    ##    ####    ##        ##       ##  ##      ##   ###
##     ##      ##   ##     ###   ##          ##   ####   #####    ###### #
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 


load("raw_data/MASTER_muATLAS_scRNAseq_cpm.RData")

####################################################################################################################################################
#FIGURE 5B:  tissue-barrier specific ABC-transporter expression from mouse cell atlas
####################################################################################################################################################


library("gplots")
library("RColorBrewer")


hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns),
            cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
            density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "log10(cpm)",key.xlab = "log10(cpm)")
}

br_idx<-unique(c(grep("endothelial",colnames(MASTER_muATLAS_cpm)),grep("enterocyte",colnames(MASTER_muATLAS_cpm)),grep("proximal",colnames(MASTER_muATLAS_cpm)),grep("hepatocyte",colnames(MASTER_muATLAS_cpm))))


hcluster(log10(t(MASTER_muATLAS_cpm[c("Abcb1a","Abcb1b","Abcg2","Abcc1","Abcc2"),br_idx]+1)),col_size = 1,row_size = 1)


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  
##   ##     ##   ##   ###     ##   ########  ##########  ##     #######   #######                                                                                                      
##   ##    ##    ##   ####    ##   ##            ##      ##    ##        ##                                                                                  
##   ##   ##     ##   ## ##   ##   ##            ##      ##   ##         ##                                                                                     
##   ######      ##   ##  ##  ##   #######       ##      ##   ##          ######                                                                                           
##   ##   ##     ##   ##   ## ##   ##            ##      ##   ##               ##                                                                               
##   ##    ##    ##   ##    ####   ##            ##      ##    ##              ##                                                                                
##   ##     ##   ##   ##     ###   #######       ##      ##     ######    ######                                                                                            
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      #######      ##      #######    #######    ##  #######   #######                                                                                    
##      ##    ##    ####     ##    ##   ##    ##   ##  ##        ##    ##                                                                                      
##      ##    ##   ##  ##    ##    ##   ##    ##   ##  ##        ##    ##                                                                                                           
##      #######   ##    ##   #######    #######    ##  #######   #######                                                                                      
##      ##    ##  ########   ##   ##    ##   ##    ##  ##        ##   ##                                                                                                       
##      ##    ##  ##    ##   ##    ##   ##    ##   ##  ##        ##    ##                                                                                                         
##      #######   ##    ##   ##     ##  ##     ##  ##  #######   ##     ##                                                                                                                        
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##        ##          ##    #######    ########     #########   ##         #########                                                              
##        ###        ###   ##     ##   ##     ##    ##          ##        ##                                                  
##        ####      ####   ##     ##   ##      ##   ##          ##        ##                                                  
##        ## ##    ## ##   ##     ##   ##      ##   #######     ##         #######                                                       
##        ##  ##  ##  ##   ##     ##   ##      ##   ##          ##               ##                                            
##        ##   ####   ##   ##     ##   ##     ##    ##          ##               ##                                            
##        ##    ##    ##    #######    ########     #########   ########   #######                                                           
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  


load("BBB_matrix.RData")

linear_fit_plot<-function(x,y,pt_colors="black",xlab=NA,ylab=NA,title=NA,xlim=NA,ylim=NA){
  #linear regression:
  lin_reg<-lm(y~x)
  #make fit equation + r-squared:
  pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
  fit_coeff<-round(coef(lin_reg),6)
  r2 <- round(summary(lin_reg)$r.squared, 2)
  rmse <- round(sqrt(mean(resid(lin_reg)^2)), 2)
  eq_r2<-paste("r = ", r2," "," "," ","p-value = ",signif(pvalue_corr,2))
  #plot data:
  plot(x,y, pch = 16, cex = 0.8, col = pt_colors,xlab=xlab,ylab=ylab,main=title,xlim=xlim,ylim=ylim)
  #add fit-line and fit-equation:
  abline(lin_reg,col=rgb(1,0,0,0.5))
  mtext(eq_r2, 3, line=-2,cex=1.5,col=rgb(1,0,0,0.5))
}

##########################################################
#CALCULATE AND SAVE GENE-FUNCTION MATRIX:  i.e. substratge specificity
##########################################################
tko_v_g2<-sort(BBB_matrix[,"TKO"]/BBB_matrix[,"abcg2 -/-"],decreasing=T)# 79 data points
tko_v_b1<-sort(BBB_matrix[,"TKO"]/BBB_matrix[,"abcb1a -/-"],decreasing=T)#79 data points

tko_v_wt<-sort(BBB_matrix[,"TKO"]/BBB_matrix[,"WT"],decreasing=T)#100 data points

name_overlap<-intersect(names(tko_v_wt),intersect(names(tko_v_b1),names(tko_v_g2)))

Function_matrix<-cbind(tko_v_g2[name_overlap],tko_v_b1[name_overlap],tko_v_wt[name_overlap])
colnames(Function_matrix)<-c("Abcb1a","Abcg2","TKO")
save(Function_matrix,file="Function_matrix.RData")


##########################################################
#FIGURE 5D:  Abcb1 vs Abcg2 specificity scatter plot + quadrants
##########################################################

linear_fit_plot(log10(tko_v_g2[name_overlap]),log10(tko_v_b1[name_overlap]),xlim=c(-0.5,2.5),ylim=c(-0.5,2.5),
                xlab="Abcb1 (tko/g2)",ylab="Abcg2 (tko/b1)",title="Abcb1-ko vs Abcg2-ko compare")
text(log10(tko_v_g2[name_overlap])+0.2,log10(tko_v_b1[name_overlap]),name_overlap,cex=0.5)
#abline(a=0,b=1,lty=2)
abline(v=mean(log10(tko_v_g2[name_overlap])),lty=2,lwd=0.5)
abline(h=mean(log10(tko_v_b1[name_overlap])),lty=2,lwd=0.5)





####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      ########   ##     ##   ########    ###      ###                                                          
##     ##          ##     ##   ##          ####    #### 
##    ##           ##     ##   ##          ## ##  ## ##             
##    ##           #########   ########    ##  ####  ##        
##    ##           ##     ##   ##          ##   ##   ##      
##     ##          ##     ##   ##          ##        ## 
##      ########   ##     ##   ########    ##        ##                    
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########    #######     ########      #######    ##              ##      ########### ########    ########   ####    ##                                                                            
##     ##          ##     ##    ##     ##     ##    ##   ##             ####         ##         ##      ##      ##  ## ##   ##                
##    ##           ##     ##    ##     ##     ##    ##   ##            ##  ##        ##         ##      ##      ##  ##  ##  ##                         
##    ##           ##     ##    ########      #######    ##           ##    ##       ##         ##      ##      ##  ##   ## ##                 
##    ##           ##     ##    ##    ##      ##   ##    ##          ##########      ##         ##      ##      ##  ##    ####                    
##     ##          ##     ##    ##     ##     ##    ##   ##         ##        ##     ##         ##      ##      ##  ##     ###                  
##      ########    #######     ##      ##    ##     ##  ########  ##          ##    ##      ########    ########   ##      ##                             
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################

#LOAD FUNCTION MATRIX (made previous section)
load("Function_matrix.RData")

#LOAD DRUGBANK DATA FOR DRUG LIBRARY
load("raw_data/DrugBank_analysis.RData")

library("gplots")
library("RColorBrewer")


hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns),
            cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
            density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "log10(cpm)",key.xlab = "log10(cpm)")
}

##########################################################
#RECONCILE DRUG BANK AND BR-data names
##########################################################
#extract b1 vector
B1_vector<-Function_matrix[,"Abcb1a"]
G2_vector<-Function_matrix[,"Abcg2"]

#renames to drugbank
names(B1_vector)<-sapply(strsplit(names(B1_vector),split="_"),function(x) x[1])
names(B1_vector)<-gsub("Abermaciclib","Abemaciclib",gsub("Axtinib","Axitinib",gsub("Flavopiridol","Alvocidib",gsub("Olamzapine","Olanzapine",gsub("Cyclosporine","Cyclosporin A",names(B1_vector))))))
names(G2_vector)<-sapply(strsplit(names(G2_vector),split="_"),function(x) x[1])
names(G2_vector)<-gsub("Abermaciclib","Abemaciclib",gsub("Axtinib","Axitinib",gsub("Flavopiridol","Alvocidib",gsub("Olamzapine","Olanzapine",gsub("Cyclosporine","Cyclosporin A",names(G2_vector))))))



#id overlap
drug_overlap<-names(B1_vector)[which(names(B1_vector) %in% rownames(DrugBank_slice))]

############################################################################################
#Figure 5D  PLOT GLOBAL CORRELATIONS
############################################################################################
linreg_data<-data.frame(cbind(log10(B1_vector[drug_overlap]),
                              log10(G2_vector[drug_overlap]),
                              as.numeric(DrugBank_slice[drug_overlap,18]),
                              as.numeric(DrugBank_slice[drug_overlap,22]),
                              as.numeric(DrugBank_slice[drug_overlap,23]),
                              as.numeric(DrugBank_slice[drug_overlap,24]),
                              as.numeric(DrugBank_slice[drug_overlap,25]),
                              as.numeric(DrugBank_slice[drug_overlap,26]),
                              as.numeric(DrugBank_slice[drug_overlap,27]),
                              as.numeric(DrugBank_slice[drug_overlap,28]),
                              as.numeric(DrugBank_slice[drug_overlap,29]),
                              as.numeric(DrugBank_slice[drug_overlap,30]),
                              as.numeric(DrugBank_slice[drug_overlap,31])))

colnames(linreg_data)<-c("Abcb1a","Abcg2","LogP","MW","MW2","Rings","Charge","pKaAcid","pKaBase","PSA","Polariz","Refract","Rotate")

sort(cor(data.matrix(linreg_data))[,"Abcb1a"],decreasing=T)
#   Abcb1a     Rotate     Charge       LogP      Abcg2    Polariz    Refract         MW        MW2      Rings    pKaBase        PSA 
#1.0000000  0.4340123  0.4052108  0.3848664  0.3554385  0.3429521  0.3420746  0.2919675  0.2914890  0.2127239  0.2114959 -0.1163568 

hcluster(t(scale(data.matrix(linreg_data),center=T,scale=T)),col_size = 0.5,title="ABC-flux vs chemical prop")


hcluster(cor(data.matrix(linreg_data),use="pairwise.complete.obs"))




