
########################
#LIBRARIES
########################
library("shiny")
library('raster')
library("gplots")
library("RColorBrewer")
library("Hmisc")

########################
#FUNCTIONS:
########################
color.bar <- function(lut, min, max=-min, title='') {
  #nticks=(round(max,0)-round(min,0)+1)/2
  nticks=5
  scale = (length(lut)-1)/(max-min)
  #ticks=seq(min, max, len=nticks)
  ticks=c(-2,-1,0,1,2)
  labels=c("0.01","0.1","1","10","100")
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=title, main='',cex.lab=1.5)
  axis(2, ticks,labels=labels, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

plot_KO_ratios<-function(drug="Doxorubicin_24hr",gene="abcb1a -/-",data){
  #drug="Doxorubicin_24hr"
  #gene="abcb1a -/-"
  
  row_idx<-intersect(which(drug_names_ordered==drug),which(gene_list_ordered==gene))

  #ADDED THIS TO MAKE KEY (11-28)
  mouse_models<-unique(gene_list_ordered[which(drug_names_ordered==drug)])
  #ADDED THIS AS SOME DRUGS HAVE MORE THAN ONE DOSE (e.g. lapatinib) (11-28)
  if(length(row_idx)>1){row_idx<-row_idx[1]}
  
  
  
  if(length(row_idx)==1){
    
    ########################
    #SET VALUES AND COLOR KEY
    ########################

    #PLUG INTO "BR vector" from original function so can leave the visualization code the same:
    BR_vector<-data[row_idx,]
    
    
    BR_vector[which(BR_vector<(0.01))]<-0.01
    BR_vector[which(BR_vector>(100))]<-100
    
    if(length(which(is.na(BR_vector)==T))>0){
      BR_vector_na<-BR_vector[-which(is.na(BR_vector)==T)]
    }else{
      BR_vector_na<-BR_vector
    }
    
    
    
    rbPal <- colorRampPalette(c('blue','white','red'))
    color_scheme_avg <- rbPal(100)[as.numeric(cut(c(log10(BR_vector_na),-2,2),breaks = 100))]
    names(color_scheme_avg)<-c(names(BR_vector_na),"min","max")
    
    
    ########################
    #RECOLOR TISSUES
    ########################
    blank_tmp<-blank_raster
    
    if ("stomach" %in% names(color_scheme_avg)){
      blank_tmp[stomach_col_idx] <- color_scheme_avg["stomach"]
      blank_tmp[stomach_blk_idx] <- 'black'
    }
    
    if ("lungs" %in% names(color_scheme_avg)){
      blank_tmp[lungs_col_idx] <- color_scheme_avg["lungs"]
      blank_tmp[lungs_blk_idx] <- 'black'
    }
    
    
    if ("heart" %in% names(color_scheme_avg)){
      blank_tmp[heart_col_idx] <- color_scheme_avg["heart"]
      blank_tmp[heart_blk_idx] <- 'black'
    }
    
    
    if ("liver" %in% names(color_scheme_avg)){
      blank_tmp[liver_col_idx] <- color_scheme_avg["liver"]
      blank_tmp[liver_blk_idx] <- 'black'
    }
    
    
    if ("intestines" %in% names(color_scheme_avg)){
      blank_tmp[intestines_col_idx] <- color_scheme_avg["intestines"]
      blank_tmp[intestines_blk_idx] <- 'black'
    }
    
    
    if ("kidney" %in% names(color_scheme_avg)){
      blank_tmp[kidney_col_idx] <- color_scheme_avg["kidney"]
      blank_tmp[kidney_blk_idx] <- 'black'
    }
    
    
    if ("brain" %in% names(color_scheme_avg)){
      blank_tmp[brain_col_idx] <- color_scheme_avg["brain"]
      blank_tmp[brain_blk_idx] <- 'black'
    }
    
    
    if ("bone" %in% names(color_scheme_avg)){
      blank_tmp[bone_col_idx] <- color_scheme_avg["bone"]
      blank_tmp[bone_blk_idx] <- 'black'
    }
    
    
    if ("muscle" %in% names(color_scheme_avg)){
      blank_tmp[muscle_col_idx] <- color_scheme_avg["muscle"]
      blank_tmp[muscle_blk_idx] <- 'black'
    }
    
    
    
    if ("testes" %in% names(color_scheme_avg)){
      blank_tmp[testes_col_idx] <- color_scheme_avg["testes"]
      blank_tmp[testes_blk_idx] <- 'black'
    }
    
    
    if ("gall" %in% names(color_scheme_avg)){
      blank_tmp[gall_col_idx] <- color_scheme_avg["gall"]
      blank_tmp[gall_blk_idx] <- 'black'
    }
    
    
    if ("spleen" %in% names(color_scheme_avg)){
      blank_tmp[spleen_col_idx] <- color_scheme_avg["spleen"]
      blank_tmp[spleen_blk_idx] <- 'black'
    }
    
    
    if ("thymus" %in% names(color_scheme_avg)){
      blank_tmp[thymus_col_idx] <- color_scheme_avg["thymus"]
      blank_tmp[thymus_blk_idx] <- 'black'
    }
    
    
    if ("lymph" %in% names(color_scheme_avg)){
      blank_tmp[lymph_col_idx] <- color_scheme_avg["lymph"]
      blank_tmp[lymph_blk_idx] <- 'black'
    }
    
    
    ########################
    #PLOT
    ########################
    raster::plot(blank_tmp,main="test")
    
    #text(x=-150, y=1500,cex=2,pos=4,labels=paste0(barrier_csv$drug[row_idx]))
    
    #text(x=-150, y=1450,cex=1.5,pos=4,labels=paste0(barrier_csv$dose[row_idx],"mg/kg (",barrier_csv$admin[row_idx],")"))
    
    #text(x=-150, y=1400,cex=1.5,pos=4,labels=paste0(barrier_csv$time[row_idx]," hrs"))
    
    #text(x=-150, y=1350,cex=1.5,pos=4,labels=paste0(barrier_csv$method[row_idx]))
    

    
    legend("topleft",title=barrier_csv$drug[row_idx],legend=c(paste0(barrier_csv$dose[row_idx],"mg/kg (",barrier_csv$admin[row_idx],")"),
                                                              barrier_csv$gene[row_idx],
                                                              paste0(barrier_csv$time[row_idx]," hrs"),
                                                              paste0(barrier_csv$method[row_idx])),cex=1,bty="n")
    
    subplot(color.bar(rbPal(100),-2,title = "WT:KO Ratio"), x=-50, y=700,size=c(0.2,5))
    
    #ADD Model Key:
    tmp_col<-c("black","black","black")
    tmp_col[-which(gene_list_4list %in% mouse_models)]<-"white"
    legend("topright",legend=gene_list_4list,text.col = tmp_col,col=tmp_col,bty="n",cex=1.5,pch=16,title="Data Available:",title.col = "black")
    
    
  }else{
    raster::plot(blank_raster)
    
    #ADD Model Key:
    tmp_col<-c("black","black","black")
    tmp_col[-which(gene_list_4list %in% mouse_models)]<-"white"
    legend("topright",legend=gene_list_4list,text.col = tmp_col,col=tmp_col,bty="n",cex=1.5,pch=16,title="Data Available:",title.col = "black")
    
    
  }
}

########################
#LOAD DATA:
########################
load("BR-mapping-GEO.RData")
barrier_csv_raw<-read.csv("ko.csv",as.is=T)
barrier_csv<-barrier_csv_raw[-which(barrier_csv_raw$gene=="WT"),]

tissues<-c("stomach","lungs","heart","liver","intestines","kidney","brain","muscle","testes","gall","spleen","thymus","lymph")
data=data.matrix(barrier_csv[,tissues])

#Added this to visualize how free concentration of drug changes:
data<-(data)^-1

########################
#MENU OPTIONS
########################
drug_names_ordered<-paste(barrier_csv$drug,paste0(barrier_csv$time,"hr"),sep="_")
drug_names_4list<-sort(unique(drug_names_ordered))

gene_list_ordered<-gsub("abcb1a/1b -/-","abcb1a -/-",barrier_csv$gene)
gene_list_4list<-c("abcb1a -/-","abcg2 -/-","TKO")



ui <- fluidPage(
  
  sidebarLayout(
    sidebarPanel(width=3,
                 selectizeInput(inputId = 'search',label = 'Select Drug Below',choices = drug_names_4list,selected="Digoxin_4hr"),
                 radioButtons(inputId = 'gene',label = 'Genetics',choices = gene_list_4list,selected = "abcb1a -/-"),
                 fileInput(inputId = 'upload',label = 'Upload CSV',accept=".csv"),
                 downloadButton('downloadPlot', 'Download Plot')),
    mainPanel(plotOutput("hist",height=800),width=9)
  )
  
)

#input & output are list-like objects
server <- function(input, output) {
  output$hist <- renderPlot({
    ########################
    #UPLOAD FILE STOP:  replace old PK DB
    ########################
    if (is.null(input$upload)==F){
      file1<-input$upload
      barrier_csv<-read.csv(file1$datapath,as.is=T)
      tissues<-c("stomach","lungs","heart","liver","intestines","kidney","brain","muscle","testes","gall","spleen","thymus","lymph")
      data=data.matrix(barrier_csv[,tissues])
    }

    
    ########################
    #INPUT PARAMETERS:
    ########################
    drug=input$search
    gene=input$gene
    

    
    ########################
    #RECOLOR TISSUES
    ########################

    plot_KO_ratios(drug,gene,data)
    
    
    
  })
  
  
  
  ####################################################################################################################################################################################
  #DOWNLOADS: #1 make plot as function (basically just copy of render-plot above)   #2 make PDF for download
  ####################################################################################################################################################################################
  
  
  
  plotInput <- function(){
    ########################
    #UPLOAD FILE STOP:  replace old PK DB
    ########################
    if (is.null(input$upload)==F){
      file1<-input$upload
      barrier_csv<-read.csv(file1$datapath,as.is=T)
      tissues<-c("stomach","lungs","heart","liver","intestines","kidney","brain","muscle","testes","gall","spleen","thymus","lymph")
      data=data.matrix(barrier_csv[,tissues])
    }
    
    
    ########################
    #INPUT PARAMETERS:
    ########################
    drug=input$search
    gene=input$gene
    
    
    
    ########################
    #RECOLOR TISSUES
    ########################
    
    plot_KO_ratios(drug,gene,data)
  }
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".pdf", sep="")
    },
    content = function(file) {
      pdf(file,width=10,height=8)
      plotInput()
      dev.off()
    },
    contentType = "image/pdf") 
  
  
}

shinyApp(ui = ui, server = server)
