library("pracma")
setwd("path/to/data")
list.files()
dir.create("output")

figfam2seed <- read.csv("figfam2SEED.csv", header=T, stringsAsFactors = F)
head(figfam2seed)

orf2figfam <- read.table("Protein2FigFAM.2col", header=F, stringsAsFactors = F)
names(orf2figfam)[1:2] <- c("FigfamID", "orf_id")
head(orf2figfam)

tax_inf <- read.table("Protein2Taxonomy.txt", header=T, stringsAsFactors = F)
head(tax_inf)

class_orf<-merge(figfam2seed, orf2figfam, by.x="FIGFAMID", by.y="FigfamID", all.x=T)
##sanity check
length(unique(class_orf$FIGFAMID))
length(unique(figfam2seed$FIGFAMID)) 
##sanity check over

class_orf_tax<-merge(class_orf, tax_inf, by.x="orf_id", by.y="orf_id", all.x=T)

#######################################################################################################
#
#
#ANALYSIS to describe total transcriptome
#
#
#
#######################################################################################################
input_files<-c("Uniq_trt_null_ind_pts1_10.FulldesignVS.1_split.txt",
               "Uniq_trt_null_ind_pts2_20.FulldesignVS.1_split.txt",
               "Uniq_trt_null_ind_pts3_30.FulldesignVS.1_split.txt",
               "Uniq_trt_null_ind_pts4_40.FulldesignVS.1_split.txt")
names(input_files)<-c("day5", "day1", "day2", "day7")
names(input_files)
combined_summary <- as.data.frame(unique(figfam2seed$Level1))
names(combined_summary) <-  "Level1"
for(i in 1:length(input_files)){
  curr_sample<-(names(input_files[i]))
  print(curr_sample)
  temp<-read.table(input_files[i], header=T, stringsAsFactors = F)
  numcol<-ncol(temp)
  temp<-temp[,(-1*(numcol-2)):(-1*numcol)]
  ncontrol<-ncol(temp[,3:8])
  ncarb<-ncol(temp[,9:(ncol(temp))])
  temp["control_sum"]<-rowSums(temp[,3:8], na.rm=T)
  temp["carbadox_sum"]<-rowSums(temp[,9:(ncol(temp))], na.rm=T)
  temp["control_avg"]<-temp$control_sum/ncontrol
  temp["carbadox_avg"]<-temp$carbadox_sum/ncarb
  temp1<-unique(temp[,3:14])
  #day2qs_names<-rownames(day2qs)
  #day2qs_names[1:5]
  print(paste("classifying", curr_sample))
  temp_SEED<-merge(class_orf_tax, temp, by.y="FigFam", by.x="FIGFAMID", all.y=T)
  temp_SEED[,3:6][is.na(temp_SEED[,3:6])]<-"Uncategorized"
  print(paste("done classifying", curr_sample))
  
  temp_summary<-matrix(ncol=5)
  line1<-c("Level1",paste(curr_sample, "control_avg_total_reads", sep=""),paste(curr_sample, "control_stderr_total_reads", sep=""),paste(curr_sample, "carb_avg_total_reads", sep=""),paste(curr_sample, "carb_stderr_total_reads", sep=""))
  line1
  colnames(temp_summary)<- line1
  #temp_summary<-t(as.matrix(temp_summary[1,1:5]))
  #temp<-as.matrix(subset(day2qs_class_q0.1_carb, Level1=="Motility___and___Chemotaxis"))
  #length(unique(temp[,1]))
  #temp[,1]
  temp_SEED$hypothetical<-regexpr('hypothetical_protein', temp_SEED$description) 
  temp_SEED[temp_SEED$hypothetical>0,]$Level1<-"Hypothetical_Protein"
  temp_SEED<-temp_SEED[,(-1*ncol(temp_SEED))]
  #print(paste("writing file for", curr_sample))
  #write.csv(temp_SEED, file=paste(curr_sample, "SEED_allFigfams.csv", sep=""))
  print(paste("counting SEED subsystems for", curr_sample))
  for (j in unique(temp_SEED$Level1)){
    tempLevel<-as.data.frame(subset(temp_SEED, Level1==j))
    tempLevel_unique<-unique(tempLevel[,c(-2,-7:-12)])
    #temp_line<-c("hi", 1,5,3,2,5)
    control_sums <- colSums(tempLevel_unique[7:(6+ncontrol)])
    carb_sums <- colSums(tempLevel_unique[(7+ncontrol):(6+ncontrol+ncarb)])
    temp_line<-c(as.character(unique(tempLevel$Level1)), sum(tempLevel_unique$control_avg), std_err(control_sums), sum(tempLevel_unique$carbadox_avg), std_err(carb_sums))
    temp_summary<-rbind(temp_summary, temp_line)
    write.csv(temp_summary, file=paste("output/", curr_sample, "SEED_allFigfams_summary.csv", sep=""))
  }
  combined_summary <- merge(combined_summary, temp_summary, all = T)
}

write.csv(combined_summary, file=paste("output/SEED_allFigfams_summary.csv", sep=""))


##########################################################################################
#
#Taxonomic Inferences
#
##########################################################################################


otu_table <- read.table("uniq.spikein.norm.txt", header=T, stringsAsFactors = F)
otu_tax_table <- merge(otu_table, tax_inf, by.x = 0, by.y = "orf_id")
colnames(otu_tax_table)[1] <- "orf_id"


###################################
#Get Counts

rm(tax_summary)
start_line <- rep(1.1, 62)
tax_summary<-data.frame(t(start_line))
samples<-colnames(otu_tax_table[2:48])
tax_line1<-c("Tax_Level","Taxa","Phylum","Class","Order","Family","Genus","day5avg_control","day5avg_carbadox","day1avg_control","day1avg_carbadox","day2avg_control","day2avg_carbadox","day7avg_control","day7avg_carbadox")
tax_line2<-c(tax_line1, samples)
tax_line2
colnames(tax_summary)<-tax_line2
#tax_summary<-t(as.matrix(tax_summary[1,]))
tax_summary
for (i in unique(otu_tax_table$Phylum)){
  tempLevel<-as.data.frame(subset(otu_tax_table, Phylum==i))
  tempLevel_sums<-colSums(tempLevel[2:48])
  tempLevel_avgs<-c("Phylum",i,tempLevel[1,50],NA,NA,NA,NA,as.numeric(mean(tempLevel_sums[c(7,15,19,25,36,40)])),as.numeric(mean(tempLevel_sums[c(3,11,23,28,32,43)])),as.numeric(mean(tempLevel_sums[c(5,13,17,34,38,45)])),as.numeric(mean(tempLevel_sums[c(1,9,21,26,30)])),as.numeric(mean(tempLevel_sums[c(6,14,18,35,39,46)])),as.numeric(mean(tempLevel_sums[c(2,10,22,27,31,42)])),as.numeric(mean(tempLevel_sums[c(8,16,20,37,41,47)])),as.numeric(mean(tempLevel_sums[c(4,12,24,29,33,44)])))
  temp_line<-c(tempLevel_avgs, tempLevel_sums)
  tax_summary<-rbind(tax_summary, temp_line)
}
str(tax_summary)

for (i in unique(otu_tax_table$Class)){
  tempLevel<-as.data.frame(subset(otu_tax_table, Class==i))
  tempLevel_sums<-colSums(tempLevel[2:48])
  tempLevel_avgs<-c("Class",i,tempLevel[1,50],tempLevel[1,51],NA,NA,NA,mean(tempLevel_sums[c(7,15,19,25,36,40)]),mean(tempLevel_sums[c(3,11,23,28,32,43)]),mean(tempLevel_sums[c(5,13,17,34,38,45)]),mean(tempLevel_sums[c(1,9,21,26,30)]),mean(tempLevel_sums[c(6,14,18,35,39,46)]),mean(tempLevel_sums[c(2,10,22,27,31,42)]),mean(tempLevel_sums[c(8,16,20,37,41,47)]),mean(tempLevel_sums[c(4,12,24,29,33,44)]))
  temp_line<-c(tempLevel_avgs, tempLevel_sums)
  tax_summary<-rbind(tax_summary, temp_line)
}

for (i in unique(otu_tax_table$Order)){
  tempLevel<-as.data.frame(subset(otu_tax_table, Order==i))
  tempLevel_sums<-colSums(tempLevel[2:48])
  tempLevel_avgs<-c("Order",i,tempLevel[1,50],tempLevel[1,51],tempLevel[1,52],NA,NA,mean(tempLevel_sums[c(7,15,19,25,36,40)]),mean(tempLevel_sums[c(3,11,23,28,32,43)]),mean(tempLevel_sums[c(5,13,17,34,38,45)]),mean(tempLevel_sums[c(1,9,21,26,30)]),mean(tempLevel_sums[c(6,14,18,35,39,46)]),mean(tempLevel_sums[c(2,10,22,27,31,42)]),mean(tempLevel_sums[c(8,16,20,37,41,47)]),mean(tempLevel_sums[c(4,12,24,29,33,44)]))
  temp_line<-c(tempLevel_avgs, tempLevel_sums)
  tax_summary<-rbind(tax_summary, temp_line)
}

for (i in unique(otu_tax_table$Family)){
  tempLevel<-as.data.frame(subset(otu_tax_table, Family==i))
  tempLevel_sums<-colSums(tempLevel[2:48])
  tempLevel_avgs<-c("Family",i,tempLevel[1,50],tempLevel[1,51],tempLevel[1,52],tempLevel[1,53],NA,mean(tempLevel_sums[c(7,15,19,25,36,40)]),mean(tempLevel_sums[c(3,11,23,28,32,43)]),mean(tempLevel_sums[c(5,13,17,34,38,45)]),mean(tempLevel_sums[c(1,9,21,26,30)]),mean(tempLevel_sums[c(6,14,18,35,39,46)]),mean(tempLevel_sums[c(2,10,22,27,31,42)]),mean(tempLevel_sums[c(8,16,20,37,41,47)]),mean(tempLevel_sums[c(4,12,24,29,33,44)]))
  temp_line<-c(tempLevel_avgs, tempLevel_sums)
  tax_summary<-rbind(tax_summary, temp_line)
}

for (i in unique(otu_tax_table$Genus)){
  tempLevel<-as.data.frame(subset(otu_tax_table, Genus==i))
  tempLevel_sums<-colSums(tempLevel[2:48])
  tempLevel_avgs<-c("Genus",i,tempLevel[1,50],tempLevel[1,51],tempLevel[1,52],tempLevel[1,53],tempLevel[1,54],mean(tempLevel_sums[c(7,15,19,25,36,40)]),mean(tempLevel_sums[c(3,11,23,28,32,43)]),mean(tempLevel_sums[c(5,13,17,34,38,45)]),mean(tempLevel_sums[c(1,9,21,26,30)]),mean(tempLevel_sums[c(6,14,18,35,39,46)]),mean(tempLevel_sums[c(2,10,22,27,31,42)]),mean(tempLevel_sums[c(8,16,20,37,41,47)]),mean(tempLevel_sums[c(4,12,24,29,33,44)]))
  temp_line<-c(tempLevel_avgs, tempLevel_sums)
  tax_summary<-rbind(tax_summary, temp_line)
}

################################

tax_summary <- tax_summary[-1,]
write.csv(tax_summary, file="output/tax_inferences_summary.csv", quote = F)
row_names_test<-(paste0("X",tax_summary[,1],"-",tax_summary[,2]))
row.names(tax_summary)<-row_names_test
phylogeny<-tax_summary[,1:7]

##Look for statistical differences:

library(QuasiSeq)

str(tax_summary)
tax_summary1 <- apply(tax_summary[,16:ncol(tax_summary)], 2, as.integer)
row.names(tax_summary1) <- row.names(tax_summary)

str(tax_summary1)
dim(tax_summary1)

day5cntl_carb<-tax_summary1[,c(7,15,19,25,36,40,3,11,23,28,32,43)]
day1cntl_carb<-tax_summary1[,c(5,13,17,34,38,45,1,9,21,26,30)]
day2cntl_carb<-tax_summary1[,c(6,14,18,35,39,46,2,10,22,27,31,42)]
day7cntl_carb<-tax_summary1[,c(8,16,20,37,41,47,4,12,24,29,33,44)]
day5cntl_carb<-day5cntl_carb[which(rowSums(day5cntl_carb)>0),]
day1cntl_carb<-day1cntl_carb[which(rowSums(day1cntl_carb)>0),]
day2cntl_carb<-day2cntl_carb[which(rowSums(day2cntl_carb)>0),]
day7cntl_carb<-day7cntl_carb[which(rowSums(day7cntl_carb)>0),]


#day5cntl_7cntl <- tax_summary1[,c(7,15,19,25,36,40,8,16,20,37,41,47)]
#day5cntl_7cntl <- day5cntl_7cntl[which(rowSums(day5cntl_7cntl)>0),]

day5cntl_carb[1:5,1:5]
log.offset=rep(1, 12)
log.offset

##Set up model
trt<-as.vector(c(rep(1,6), rep(2,6)))
mn<-as.vector(rep(1,12))
trt
design.list<-vector('list',2)
design.list[[1]]<-model.matrix(~trt)
design.list[[2]]<-mn
is.na(day5cntl_carb)<-0
test<-as.matrix(day5cntl_carb[1:5,1:5])
test1<-as.numeric(test)
test2<-round(as.numeric(day5cntl_carb[1:5,1:5]))
dim(test2)

##Run QuasiSeq
fit2_day5<-QL.fit(day5cntl_carb, design.list,log.offset=log.offset, Model='NegBin')
fit2_day2<-QL.fit(day2cntl_carb, design.list,log.offset=log.offset, Model='NegBin')
fit2_day7<-QL.fit(day7cntl_carb, design.list,log.offset=log.offset, Model='NegBin')
#fit2_day5_7<-QL.fit(day5cntl_7cntl, design.list,log.offset=log.offset, Model='NegBin')

trt<-as.vector(c(rep(1,6), rep(2,5)))
mn<-as.vector(rep(1,11))
log.offset=rep(1, 11)
trt
design.list<-vector('list',2)
design.list[[1]]<-model.matrix(~trt)
design.list[[2]]<-mn
fit2_day1<-QL.fit(day1cntl_carb, design.list,log.offset=log.offset, Model='NegBin')

# fit2<-QL.fit(t(phage_2000_subsample2_D2D7_gene), design.list,log.offset=log.offset, Model='NegBin')
# fit<-QL.fit(t(phage_2000_subsample2_D2D7_gene), design.list,log.offset=log.offset, Model='NegBin')
# warnings()

##Get Results
#Day 5 time-matched
results<-QL.results(fit2_day5)
results$Q.values$QLSpline
apply(results$Q.values[[3]]<0.1,2,sum)
apply(results$P.values[[3]]<0.01,2,sum)

qs  <-results$Q.values[[3]]
ps <-results$P.values[[3]]
ps_qs<-merge(ps, qs, by=0)
colnames(ps_qs)<-c("taxa","pvalue", "qvalue")
day5cntl_carb_ps_qs<-merge(day5cntl_carb, ps_qs, by.x=0, by.y=1)
day5cntl_carb_ps_qs1<-merge(day5cntl_carb_ps_qs, phylogeny, by.x=1, by.y=0, all.x = T)
write.csv(day5cntl_carb_ps_qs1, file="output/day5cntl_carb_ps_qs1.csv", quote=F)

#Day 1 time-matched
results<-QL.results(fit2_day1)
results$Q.values$QLSpline
apply(results$Q.values[[3]]<0.1,2,sum)
apply(results$P.values[[3]]<0.01,2,sum)

qs  <-results$Q.values[[3]]
ps <-results$P.values[[3]]
ps_qs<-merge(ps, qs, by=0)
colnames(ps_qs)<-c("taxa","pvalue", "qvalue")
day1cntl_carb_ps_qs<-merge(day1cntl_carb, ps_qs, by.x=0, by.y=1)
day1cntl_carb_ps_qs1<-merge(day1cntl_carb_ps_qs, phylogeny, by.x=1, by.y=0, all.x = T)
write.csv(day1cntl_carb_ps_qs1, file="output/day1cntl_carb_ps_qs1.csv", quote=F)

#Day 2 time-matched
results<-QL.results(fit2_day2)
results$Q.values$QLSpline
apply(results$Q.values[[3]]<0.1,2,sum)
apply(results$P.values[[3]]<0.01,2,sum)

qs  <-results$Q.values[[3]]
ps <-results$P.values[[3]]
ps_qs<-merge(ps, qs, by=0)
colnames(ps_qs)<-c("taxa","pvalue", "qvalue")
day2cntl_carb_ps_qs<-merge(day2cntl_carb, ps_qs, by.x=0, by.y=1)
day2cntl_carb_ps_qs1<-merge(day2cntl_carb_ps_qs, phylogeny, by.x=1, by.y=0, all.x = T)
write.csv(day2cntl_carb_ps_qs1, file="output/day2cntl_carb_ps_qs1.csv", quote=F)

#Day 7 time-matched
results<-QL.results(fit2_day7)
results$Q.values$QLSpline
apply(results$Q.values[[3]]<0.1,2,sum)
apply(results$P.values[[3]]<0.02,2,sum)

qs  <-results$Q.values[[3]]
ps <-results$P.values[[3]]
ps_qs<-merge(ps, qs, by=0)
colnames(ps_qs)<-c("taxa","pvalue", "qvalue")
day7cntl_carb_ps_qs<-merge(day7cntl_carb, ps_qs, by.x=0, by.y=1)
day7cntl_carb_ps_qs1<-merge(day7cntl_carb_ps_qs, phylogeny, by.x=1, by.y=0, all.x = T)
write.csv(day7cntl_carb_ps_qs1, file="output/day7cntl_carb_ps_qs1.csv", quote=F)

#results<-QL.results(fit2_day5_7)
#results$Q.values$QLSpline
#apply(results$Q.values[[3]]<0.1,2,sum)
#apply(results$P.values[[3]]<0.02,2,sum)

#qs  <-results$Q.values[[3]]
#ps <-results$P.values[[3]]
#ps_qs<-merge(ps, qs, by=0)
#colnames(ps_qs)<-c("taxa","pvalue", "qvalue")
#day5cntl_7cntl_ps_qs<-merge(day5cntl_7cntl, ps_qs, by.x=0, by.y=1)
#day5cntl_7cntl_ps_qs1<-merge(day5cntl_7cntl_ps_qs, phylogeny, by.x=1, by.y=0, all.x = T)
#day5cntl_7cntl_ps_qs1_0.1 <- as.data.frame(subset(day5cntl_7cntl_ps_qs1, qvalue<0.1))
#write.csv(day7cntl_carb_ps_qs1, file="day7cntl_carb_ps_qs1.csv", quote=F)



###############################################

#Analysis of Day 2 DE FigFams

###############################################



day2qs<-read.table("Uniq_trt_null_ind_pts3_30.FulldesignVS.1_split.txt", header=T)
day2qs["control"]<-rowSums(day2qs[,3:8], na.rm=T)
day2qs["Carbadox"]<-rowSums(day2qs[,9:14], na.rm=T)
day2qs["EnrichedIn"]<-ifelse(day2qs$control>day2qs$Carbadox, "Control", "Carbadox")


day2qs_SEED<-merge(class_orf_tax, day2qs, by.y="FigFam", by.x="FIGFAMID", all=T)
day2qs_SEED[,3:6][is.na(day2qs_SEED[,3:6])]<-"Uncategorized"
day2qs_class_q0.1<-as.data.frame(subset(day2qs_SEED, Qvalues<0.1))
day2qs_class_q0.1_control<-as.data.frame(subset(day2qs_class_q0.1, EnrichedIn=="Control"))
day2qs_class_q0.1_carb<-as.data.frame(subset(day2qs_class_q0.1, EnrichedIn=="Carbadox"))
carb_summary<-matrix(ncol=4)
line1<-c("Level1","carb_enriched_count","carb_avg_FC","carb_std_err_FC")
colnames(carb_summary) <-line1

day2qs_class_q0.1_carb$hypothetical<-regexpr('hypothetical_protein', day2qs_class_q0.1_carb$description) 
day2qs_class_q0.1_carb[day2qs_class_q0.1_carb$hypothetical>0,]$Level1<-"Hypothetical_Protein"
min(day2qs_class_q0.1_carb$fold_change)
max(day2qs_class_q0.1_control$fold_change)
for (i in unique(day2qs_class_q0.1_carb$Level1)){
  temp<-as.data.frame(subset(day2qs_class_q0.1_carb, Level1==i))
  #temp_line<-c("hi", 1,5,3,2,5)
  temp_line<-c(as.character(unique(temp$Level1)), length(unique(temp[,1])), mean(unique(temp$fold_change)), std_err(unique(temp$fold_change)))
  carb_summary<-rbind(carb_summary, temp_line)
}




control_summary<-matrix(ncol=4)
line1<-c("Level1","carb_decreased_count","carb_avg_FC","carb_std_err_FC")
colnames(control_summary) <- line1

day2qs_class_q0.1_control$hypothetical<-regexpr('hypothetical_protein', day2qs_class_q0.1_control$description) 
day2qs_class_q0.1_control[day2qs_class_q0.1_control$hypothetical>0,]$Level1<-"Hypothetical_Protein"
for (i in unique(day2qs_class_q0.1_control$Level1)){
  temp<-as.data.frame(subset(day2qs_class_q0.1_control, Level1==i))
  #temp_line<-c("hi", 1,5,3,2,5)
  temp_line<-c(as.character(unique(temp$Level1)), length(unique(temp[,1])), mean(unique(temp$fold_change)), std_err(unique(temp$fold_change))) 
  control_summary<-rbind(control_summary, temp_line)
}

carb_control_summary<-merge(carb_summary,control_summary, by="Level1", all=T)

write.csv(carb_control_summary, file="output/Carb_Control_enriched_Figfams.csv")
write.csv(day2qs_class_q0.1, file="output/day2qs_classified_q0.1.csv")




######################################################################################################
#Day 7 analysis
#
#
#######################################################################################################

day7qs<-read.table("Uniq_trt_null_ind_pts4_40.FulldesignVS.1_split.txt", header=T)
day7qs["control"]<-rowSums(day7qs[,3:8], na.rm=T)
day7qs["Carbadox"]<-rowSums(day7qs[,9:14], na.rm=T)
day7qs["EnrichedIn"]<-ifelse(day7qs$control>day7qs$Carbadox, "Control", "Carbadox")

day7qs_SEED<-merge(class_orf_tax, day7qs, by.y="FigFam", by.x="FIGFAMID", all.y=T)
day7qs_SEED[,3:6][is.na(day7qs_SEED[,3:6])]<-"Uncategorized"
day7qs_class_q0.1<-as.data.frame(subset(day7qs_SEED, Qvalues<0.1))
day7qs_class_q0.1_control<-as.data.frame(subset(day7qs_class_q0.1, EnrichedIn=="Control"))
day7qs_class_q0.1_carb<-as.data.frame(subset(day7qs_class_q0.1, EnrichedIn=="Carbadox"))
carb_summary7<-matrix(ncol=4)
line1<-c("Level1","carb_enriched_count","carb_avg_FC","carb_std_err_FC")
colnames(carb_summary7) <- line1
carb_summary7<-t(as.matrix(carb_summary7[1,1:4]))

day7qs_class_q0.1_carb$hypothetical<-regexpr('_hypothetical_protein_', day7qs_class_q0.1_carb$description) 
day7qs_class_q0.1_carb[day7qs_class_q0.1_carb$hypothetical>0,]$Level1<-"Hypothetical_Protein"
min(day7qs_class_q0.1_carb$fold_change)
max(day7qs_class_q0.1_control$fold_change)
for (i in unique(day7qs_class_q0.1_carb$Level1)){
  temp<-as.data.frame(subset(day7qs_class_q0.1_carb, Level1==i))
  #temp_line<-c("hi", 1,5,3,2,5)
  temp_line<-c(as.character(unique(temp$Level1)), length(unique(temp[,1])), mean(unique(temp$fold_change)), std_err(unique(temp$fold_change))) 
  carb_summary7<-rbind(carb_summary7, temp_line)
}



control_summary7<-matrix(ncol=4)
line1<-c("Level1","carb_decreased_count","carb_avg_FC","carb_std_err_FC")
colnames(control_summary7) <- line1
control_summary7<-t(as.matrix(control_summary7[1,1:4]))

day7qs_class_q0.1_control$hypothetical<-regexpr('_hypothetical_protein_', day7qs_class_q0.1_control$description) 
day7qs_class_q0.1_control[day7qs_class_q0.1_control$hypothetical>0,]$Level1<-"Hypothetical_Protein"
for (i in unique(day7qs_class_q0.1_control$Level1)){
  temp<-as.data.frame(subset(day7qs_class_q0.1_control, Level1==i))
  #temp_line<-c("hi", 1,5,3,2,5)
  temp_line<-c(as.character(unique(temp$Level1)), length(unique(temp[,1])), mean(unique(temp$fold_change)), std_err(unique(temp$fold_change))) 
  control_summary7<-rbind(control_summary7, temp_line)
}

carb_control_summary7<-merge(carb_summary7,control_summary7, by="Level1", all=T)

write.csv(carb_control_summary7, file="output/Carb_Control_enriched_Figfams_stats_day7_20160502.csv")
write.csv(day7qs_class_q0.1, file="output/day7qs_classified_q0.1_20160502")
write.csv(day7qs_class_q0.1_carb, file="output/day7qs_carb_classified_q0.1_20160502")
write.csv(day7qs_class_q0.1_control, file="output/day7qs_control_classified_q0.1_20160502")


#######################################################################################################
library(ggplot2)
library(vegan)
library(MASS)

#Getting into shape for ordination plots
day5qs<-read.table("Uniq_trt_null_ind_pts1_10.FulldesignVS.1_split.txt", header=T)
day1qs<-read.table("Uniq_trt_null_ind_pts2_20.FulldesignVS.1_split.txt", header=T)
day2qs<-read.table("Uniq_trt_null_ind_pts3_30.FulldesignVS.1_split.txt", header=T)
day7qs<-read.table("Uniq_trt_null_ind_pts4_40.FulldesignVS.1_split.txt", header=T)



day2qs_ord<-t(day2qs[3:14])
day7qs_ord<-t(day7qs[3:14])
day5qs_ord<-t(day5qs[3:14])
day1qs_ord<-t(day1qs[3:13])
day5_1<-merge(day5qs[1:14], day1qs[1:13], all=T, by.x=c(0:2), by.y=(0:2))
day5_1_2<-merge(day5_1, day2qs[1:14], all=T, by.x=c(1:3), by.y=c(0:2))
day5_1_2_7<-merge(day5_1_2, day7qs[1:14], all=T, by.x=c(1:3), by.y=c(0:2))
day5_1_2_7[is.na(day5_1_2_7)]<-0
row.names(day5_1_2_7)<-day5_1_2_7[,1]
day5_1_2_7<-day5_1_2_7[,-1]
#write.csv(day5_1_2_7, file = "diversity_analysis_mothur/day5_1_2_7.csv", row.names = T)

day5_1_2_7_SEED<-merge(figfam2seed, day5_1_2_7, by.y="FigFam", by.x="FIGFAMID", all.y=T)
day5_1_2_7_SEED[,2:5][is.na(day5_1_2_7_SEED[,2:5])]<-"Uncategorized"
day5_1_2_7_SEED_uncategorized<-as.data.frame(subset(day5_1_2_7_SEED, Level1=="Uncategorized"))

carbadox <- c(rep(15, 6), rep(20, 6), rep(15, 6), rep(20, 5), rep(15, 6), rep(20, 6), rep(15, 6), rep(20, 6))
trt<-c(rep("nonD5",6), rep("carbD5",6), rep("nonD1",6), rep("carbD1",5), rep("nonD2",6), rep("carbD2",6), rep("nonD7",6), rep("carbD7",6))
trt
mycolors <- c(rep('#e41a1c', 12), rep('#377eb8', 11), rep('#4daf4a', 12), rep('#984ea3', 12))
metaT.mds <- metaMDS(t(day5_1_2_7[,3:ncol(day5_1_2_7)]), trace = T, display ="sites", autotransform = T)
metaT.mds
#plot.new()
#head(metaT.mds$points)
metaT.mds.points<-as.data.frame(metaT.mds$points)
ggplot(metaT.mds.points, aes(MDS1, MDS2)) +
  geom_point(colour=mycolors, shape = carbadox, show.legend=T, cex=4)
ordiplot(metaT.mds, type="p", display = "sites")
for (i in unique(trt)){
  ordihull(metaT.mds$point[grep(i,trt),], draw= "polygon", groups=trt[trt==i], col =mycolors[grep(i,trt)], label=T)
}

ord.fil_all <- envfit(metaT.mds~carbadox, perm=999)

alldaysnmds<- as.data.frame(metaT.mds$points)




#######################################################################################################
#Getting into shape for ordination plots
library(stats)
library(vegan)
trt<-c(1,1,1,1,1,1,2,2,2,2,2,2)
day2qs<-read.table("Uniq_trt_null_ind_pts3_30.FulldesignVS.1_split.txt", header=T)
day2qs_ord<-t(day2qs[3:14])
dim(day2qs_ord)
length(day2qs$FigFam)
length(unique(day2qs$FigFam))
day2qs_ord.env <- as.data.frame(row.names(day2qs_ord))
row.names(day2qs_ord.env) <- row.names(day2qs_ord)
day2qs_ord.env$treatment <- c(rep("non", 6), rep("carb", 6))
attach(day2qs_ord.env)


metaT.mds <- metaMDS(day2qs_ord, trace = T, display ="sites", autotransform = T)
metaT.mds
metaT.mds$points
plot.new()
plot(metaT.mds, type="n", display = "sites")
#ordispider(metaT.mds, treatment, col = "blue", label = TRUE)
ordiellipse(metaT.mds, treatment, conf = .95, kind = 'se', col = "black", lwd = 2)
points(metaT.mds, pch=16, col= c(rep("blue", 6), rep("red", 6)), cex=1.3)
#DON'T DO THE FOLLOWING LINE. It puts all the 12000 species on the plot.
#text(metaT.mds, display = "species", cex = 0.7, col = "black")

ord.fit <- envfit(metaT.mds~treatment, data = day2qs_ord.env, perm = 1000000)
ord.fit
plot(ord.fit)


#######################################################################################################


