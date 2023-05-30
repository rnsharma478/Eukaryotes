
library(data.table)

tss_df <- read.csv("D:\\kopal_iitd\\Eukaryotes\\Trimer\\combined_TSS_eukaryotes_tri.csv",header = TRUE, stringsAsFactors = FALSE)
#tss_transpose <-tss_df
rownames(tss_df) <- tss_df$X
tss_df$X <- NULL
tss_transpose <- transpose(tss_df)
rownames(tss_transpose) <- colnames(tss_df)
colnames(tss_transpose) <- rownames(tss_df)

cds_df <- read.csv("D:\\kopal_iitd\\Eukaryotes\\Trimer\\combine_CDS_eukaryotes_tri.csv",header = TRUE, stringsAsFactors = FALSE)
#cds_transpose <- cds_df
rownames(cds_df) <- cds_df$X
cds_df$X <- NULL
cds_transpose <- transpose(cds_df)
rownames(cds_transpose) <- colnames(cds_df)
colnames(cds_transpose) <- rownames(cds_df)


#the parameter list a,b,c,d... is not with mapping (synced) dinucleotide. It is directly in order. The reference can be found in "param ref for graph plotting" jpeg. Later for training the paramters were encoded according to dinucleotide 
lst_paras <- list('a','b','c','d','f','g','h','i','j','k','t','u','v','w','x','y','z','aa','ab','ac','ad','ae')
lst_paras_cds <- list('a','b','c','d','f','g','h','i','j','k','t','u','v','w','x','y','z','aa','ab','ac','ad','ae')

for (z in 1:length(lst_paras)){
  v <- c(unlist(tss_transpose[lst_paras[[z]]]))
  t <- c(unlist(cds_transpose[lst_paras_cds[[z]]]))
  
  Mint <- min(t)
  Maxt <- max(t)
  Minv <- min(v)
  Maxv <- max(v)
  
  minimum <- min(Mint, Maxt, Minv, Maxv)
  maximum <- max(Mint, Maxt, Minv, Maxv)
  
  
  lst_paras_expanded <- list("X Displacement","Y Displacement", "Inclination","Tip", "Shear", "Stretch", "Stagger", "Buckle", "Propel", "Opening","Alpha","Beta","Gamma","Delta","Epsil","Zeta","Chi", "Phase","Ampli","Hydrogen Bond","Stacking Energy","Solvation")
  
  # Plot the bar chart.
  
  png(sprintf("D:\\kopal_iitd\\Eukaryotes\\Trimer\\graphs of eukaryotes trimer\\%s.png",lst_paras_expanded[z]), width = 700 , height = 600)
  plot(x = -500:474,t,type = 'l', col = "blue", xlab = "Sequence Length", ylab = sprintf("%s",lst_paras_expanded[z]) , main = "combined trimer", ylim= c(minimum,maximum))
  lines(x = -500:474,v, type = "l", col = "red")
  grid(nx = NULL, ny = 8, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
  legend("bottomright", legend=c("TSS", "CDS"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
  
  dev.off()
}

