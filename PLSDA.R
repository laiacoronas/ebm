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
metadata<-dplyr::select(file, c("RecordID","SubjectID","Group",
                                           "Age","Gender","Side",
                                           "DataSet","Filename"  ,  "RaterA_Video"  ,
                                           "RaterB_Video","RaterC_Video","RaterD_Video",
                                            "RaterE_Video","RaterF_Video","RaterA_3DAvatar",
                                            "RaterB_3DAvatar","RaterC_3DAvatar","RaterD_3DAvatar",
                                           "RaterE_3DAvatar","RaterF_3DAvatar","Classifier_J48" ))
data_for_plsda<-select(file, -c("RecordID","SubjectID","Group",
                                                "Age","Gender","Side",
                                                "DataSet","Filename"  ,  "RaterA_Video"  ,
                                                "RaterB_Video","RaterC_Video","RaterD_Video",
                                                "RaterE_Video","RaterF_Video","RaterA_3DAvatar",
                                                "RaterB_3DAvatar","RaterC_3DAvatar","RaterD_3DAvatar",
                                                "RaterE_3DAvatar","RaterF_3DAvatar","Classifier_J48" ))
#Labels for the PLSDA
metadata$Labels<-rep(0,nrow(metadata))
metadata$Labels[metadata$Group=="PD patient"]<-1
metadata$Labels[metadata$Group=="Control"]<-0

#Split Training and Test, please make a better split
Train_data<-data_for_plsda[1:80,]
Train_metadata<-metadata[1:80,]
Test_data<-dplyr::setdiff(data_for_plsda,Train_data)
Test_metadata<-dplyr::setdiff(metadata,Train_metadata)#Scale data
Train_data<-as.matrix((Train_data))
#Selection of LV
num_pcs<-4
# Training a model
Labels_Train<-Train_metadata$Labels
Labels_Test<-Test_metadata$Labels
plsda_model <-mixOmics::plsda(Train_data,Labels_Train,ncomp = num_pcs)
scores_pca<-as.data.frame(plsda_model$variates$X)
#Predict Test
Test_data<-as.matrix(Test_data)
Test_predicted<-predict(plsda_model,Test_data)
Predicted_Labels<-unlist(Test_predicted$predict[,1,4])
#You can set the sensitivity of the model
Predicted_Labels[Predicted_Labels<=0.5]<-0
Predicted_Labels[Predicted_Labels>=0.5]<-1

