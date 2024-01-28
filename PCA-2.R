setwd("/Users/laiacoronassala/Desktop")
library("ggplot2")
library("dplyr")
library("tidyr")
library("mixOmics")
library("MASS")
library("sfsmisc")
library("e1071")
library("class")
library("caret")
library("plotly")
file<-readr::read_csv("parkinsonresearch_metadata_pronation-supination.csv",na = "?")
metadata<-select(file, c("RecordID","SubjectID","Group",
                                           "Age","Gender","Side",
                                           "DataSet","Filename"  ,  "RaterA_Video"  ,
                                           "RaterB_Video","RaterC_Video","RaterD_Video",
                                            "RaterE_Video","RaterF_Video","RaterA_3DAvatar",
                                            "RaterB_3DAvatar","RaterC_3DAvatar","RaterD_3DAvatar",
                                           "RaterE_3DAvatar","RaterF_3DAvatar","Classifier_J48" ))
data_for_pca<-select(file, -c("RecordID","SubjectID","Group",
                                                "Age","Gender","Side",
                                                "DataSet","Filename"  ,  "RaterA_Video"  ,
                                                "RaterB_Video","RaterC_Video","RaterD_Video",
                                                "RaterE_Video","RaterF_Video","RaterA_3DAvatar",
                                                "RaterB_3DAvatar","RaterC_3DAvatar","RaterD_3DAvatar",
                                                "RaterE_3DAvatar","RaterF_3DAvatar","Classifier_J48" ))
num_pcs<-4
data_for_pca<-as.matrix((data_for_pca))
data_for_pca_scaled <-scale(data_for_pca,center = T,scale=T)
pca_model <-mixOmics::pca(data_for_pca_scaled,ncomp = num_pcs)
scores_pca<-as.data.frame(pca_model$variates)

colnames(scores_pca) <- paste0("LV", 1:num_pcs)
PCA_dataframe <- cbind(metadata,
                       scores_pca)
varpercent =100*pca_model$prop_expl_var$X

ggplot(PCA_dataframe)+
  # geom_point(aes(x=PC1, y=PC2))+
  geom_point(aes(x=LV1, y=LV2,color=as.factor(Group)))+
  # geom_point(aes(x=PC1, y=PC2,color=as.factor(Behavior), shape=as.factor(Genotype)))+
  # ylim(-5.1,3)+
  # scale_color_manual(values=colors1)+
  xlab(sprintf("PC1 (%4.2f%%)", varpercent[1])) +
  ylab(sprintf("PC2 (%4.2f%%)", varpercent[2]))+ 
  #ggtitle("Experiments") + theme(axis.text=element_text(size=14),
  #                                axis.title=element_text(size=18,face="bold"),
  #                                plot.title = element_text(color="black", size=10),
  #                                legend.text=element_text(size=16))+
  theme(legend.title=element_blank())
fig <- plot_ly(PCA_dataframe, x = ~LV1, y = ~LV2, z = ~LV3, color = ~PCA_dataframe$Group) %>%
  add_markers(size = 12)


fig <- fig %>%
  layout(
    title = "tit",
    scene = list(bgcolor = "#e5ecf6")
  )

fig
