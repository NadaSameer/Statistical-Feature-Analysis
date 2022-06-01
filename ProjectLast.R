install.packages("stringr")
library("stringr")
library("dplyr")
#Q1

proteomes_data <- read.csv("77_cancer_proteomes_CPTAC_itraq.csv")
clinical_data <- read.csv("clinical_data_breast_cancer.csv")
PAM_data<- read.csv("PAM50_proteins.csv")

#1.1
PAM_data_tobejoined<-as.data.frame(PAM_data$RefSeqProteinID)
colnames(PAM_data_tobejoined)[1]<-"RefSeq_accession_number"

proteomes_data<-inner_join(proteomes_data, PAM_data_tobejoined, by = "RefSeq_accession_number")

clinical_id_record <- c()
for(i in 1:length(clinical_data$Complete.TCGA.ID))
{
  id_name <- clinical_data$Complete.TCGA.ID[i]
  modefy_id <- substr(clinical_data$Complete.TCGA.ID[i],6,nchar(id_name))
  modefy_id <- str_replace(modefy_id, '-' ,'.')
  clinical_data$Complete.TCGA.ID[i] <- str_replace(clinical_data$Complete.TCGA.ID[i], id_name, modefy_id)
  clinical_id_record[i] <- clinical_data$Complete.TCGA.ID[i]
}

proteomes_col <- c()
for(j in 4:length(colnames(proteomes_data)))
{
  p_id <- colnames(proteomes_data)[j]
  modefy_p_id <- substr(p_id,1,nchar(p_id)-7)
  colnames(proteomes_data)[j] <- str_replace(colnames(proteomes_data)[j], p_id, modefy_p_id)
  proteomes_col[j] <- colnames(proteomes_data)[j]
}
proteomes_col <- as.character(na.omit(proteomes_col))
proteomes_col

common <- intersect(proteomes_col , clinical_id_record)
common <- unique(append(colnames(proteomes_data)[1],common))

processed_data <- subset(proteomes_data, select = common)

PAM_data_tobejoined<-as.data.frame(PAM_data$RefSeqProteinID)
colnames(PAM_data_tobejoined)[1]<-"RefSeq_accession_number"

data_final<-inner_join(processed_data, PAM_data_tobejoined, by = "RefSeq_accession_number")


#1.2

final_data<- na.omit(data_final)

col_names_after_transport<-final_data$RefSeq_accession_number
transposed_data<-as.data.frame(t(final_data[,2:78]))
colnames(transposed_data)<-col_names_after_transport
transposed_data <- cbind(rownames(transposed_data), data.frame(transposed_data, row.names=NULL))
colnames(transposed_data)[1] <- "Complete.TCGA.ID"
clinical_data_subset <-clinical_data [,c("Complete.TCGA.ID" , 'HER2.Final.Status')]
protein_table <-  inner_join(clinical_data_subset, transposed_data, by = "Complete.TCGA.ID")

#Q2
protein_table_copy<-protein_table

#2.1
her_status_numeric<- unclass(as.factor(protein_table$HER2.Final.Status))
protein_table$HER2.Final.Status<-as.numeric(her_status_numeric)

correlation<-c()
correlation_copy<-c()
protein_info<-c()
for (i in 3:ncol(protein_table))
{
  cor_col<-cor(protein_table$HER2.Final.Status,as.numeric(protein_table[,i]),method = "pearson")
  correlation<-append(correlation,cor_col)
  correlation_copy<-append(correlation_copy,cor_col)
  protein_info<-append(protein_info,names(protein_table)[i])
}

#2.2
abs_correlation<-abs(correlation)
sorted_correlation<-sort(abs_correlation,decreasing = TRUE)
protein_info_sorted<-c()
actual_sorted<-c()
for(i in 1:length(correlation_copy))
{
  for(j in 1:length(sorted_correlation))
  {
    if(abs(correlation_copy[j])==sorted_correlation[i])
    {
      actual_sorted<-append(actual_sorted,correlation_copy[j])
      protein_info_sorted<-append(protein_info_sorted,protein_info[j])
    }
  }
}
abs_correlation
sorted_correlation

#2.3
par(mfrow=c(1,2))
plot(sorted_correlation)
plot(actual_sorted)
threshold<-0.1
listed_threshold<-c()
liste_threshold_protein_info<-c()
for(i in 1:length(sorted_correlation))
{
  if(sorted_correlation[i]>=threshold)
  {
    for(j in 1:length(sorted_correlation))
    {
      if(abs(correlation_copy[j])==sorted_correlation[i])
      {
        listed_threshold<-append(listed_threshold,correlation_copy[j])
        liste_threshold_protein_info<-append(liste_threshold_protein_info,protein_info[j])
      }
    }
  }
}
protein_info_sorted_threshold_data_frame<-data.frame(liste_threshold_protein_info,listed_threshold)
names(protein_info_sorted_threshold_data_frame)[1]<-"protein_name"
names(protein_info_sorted_threshold_data_frame)[2]<-"correlation"

#Q3
#3.1
positive_data<-protein_table_copy[protein_table_copy$HER2.Final.Status=="Positive",]
negative_data<-protein_table_copy[protein_table_copy$HER2.Final.Status=="Negative",]

#3.2
test<-c()
test_copy<-c()
test_info<-c()
test_p_value_accepted<-c()
test_p_value_accepted_info<-c()
test_p_value_rejected<-c()
test_p_value_rejected_info<-c()

#null hypothesis: there is no difference
#alternative hypothesis: there is a difference

for(i in 3:ncol(negative_data))
{
  t_test <-t.test(negative_data[,i],positive_data[,i])
  value<-t_test$statistic # return the t value test
  test<-append(test,value)
  test_copy<-append(test_copy,value)
  test_info<-append(test_info,names(negative_data)[i])
  if( t_test$p.value > 0.05 ) #accept null hypothesis
  {
    test_p_value_accepted<-append(test_p_value_accepted,t_test$p.value)
    test_p_value_accepted_info<-append(test_p_value_accepted_info,names(negative_data)[i])
  }
  else
  {
    test_p_value_rejected<-append(test_p_value_rejected,t_test$p.value)
    test_p_value_rejected_info<-append(test_p_value_rejected_info,names(negative_data)[i])
  }


}
test_p_value_accepted_info_data_frame<-data.frame(test_p_value_accepted_info,test_p_value_accepted)
names(test_p_value_accepted_info_data_frame)[1]<-"protein_name_accepted"
names(test_p_value_accepted_info_data_frame)[2]<-"p_value_accepted"
print(paste("the protein names that is accepted: ",test_p_value_accepted_info_data_frame$protein_name))

test_p_value_rejected_info_data_frame<-data.frame(test_p_value_rejected_info,test_p_value_rejected)
names(test_p_value_rejected_info_data_frame)[1]<-"protein_name_rejected"
names(test_p_value_rejected_info_data_frame)[2]<-"p_value_rejected"
print(paste("the protein names that is rejected: ",test_p_value_rejected_info_data_frame$protein_name))

#3.3
abs_test<-abs(test)
sorted_test<-sort(abs_test,decreasing = TRUE)
test_info_sorted<-c()
actual_test_sorted<-c()
for(i in 1:length(test_copy))
{
  for(j in 1:length(sorted_test))
  {
    if(abs(test_copy[j])==sorted_test[i])
    {
      actual_test_sorted<-append(actual_test_sorted,test_copy[j])
      test_info_sorted<-append(test_info_sorted,test_info[j])
    }
  }
}

test_info_sorted_data_frame<-data.frame(actual_test_sorted,test_info_sorted)
par(mfrow=c(1,2))
plot(sorted_test)
plot(actual_test_sorted)

threshold_test<-1.1
test_threshold<-c()
test_threshold_protein_info<-c()
for(i in 1:length(sorted_test))
{
  if(sorted_test[i]>=threshold_test)
  {
    for(j in 1:length(sorted_test))
    {
      if(abs(test_copy[j])==sorted_test[i])
      {
        test_threshold<-append(test_threshold,test_copy[j])
        test_threshold_protein_info<-append(test_threshold_protein_info,test_info[j])
      }
    }
  }
}

test_threshold_protein_info_data_frame<-data.frame(test_threshold_protein_info,test_threshold)
names(test_threshold_protein_info_data_frame)[1]<-"protein_name"
names(test_threshold_protein_info_data_frame)[2]<-"t_test"
#3.4
#if they are the same length
founded_in2<-c()
corr_founded<-c()
test_founded<-c()
for (i in 1:nrow(protein_info_sorted_threshold_data_frame))
{
  for(j in 1:nrow(test_threshold_protein_info_data_frame))
  {
    if(protein_info_sorted_threshold_data_frame[i,"protein_name"]==test_threshold_protein_info_data_frame[j,"protein_name"])
    {
      founded_in2<-append(founded_in2,protein_info_sorted_threshold_data_frame[i,"protein_name"])
      corr_founded<-c(corr_founded,protein_info_sorted_threshold_data_frame[i,"correlation"])
      test_founded<-c(test_founded,test_threshold_protein_info_data_frame[j,"t_test"])
    }
  }
}

#i took threshold that gives 19 protein ids so when i compare the protein ids of correlation and t test i found that 18 protein matched 
# and only one unmatched protein id 
#not similar result so they don't select the same features
#matched:  "NP_004439"    "NP_001159403" "NP_058519"    "NP_008950"    "NP_001116539" "NP_058518"    "NP_000917"    "NP_000116"    "NP_000413"    "NP_065178"    "NP_077006"    "NP_001035932" "NP_000415"    "NP_000517"    "NP_004487"    "NP_003003"    "NP_054895"    "NP_005931"
print(paste("the same features from correlation and t test: ",founded_in2,corr_founded,test_founded))
print(paste("the different feature from correlation and t tes: ",setdiff(test_threshold_protein_info_data_frame$protein_name,founded_in2),setdiff(protein_info_sorted_threshold_data_frame$correlation,corr_founded),setdiff(test_threshold_protein_info_data_frame$t_test,test_founded)))

