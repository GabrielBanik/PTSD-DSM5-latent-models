dataCFA <- data

table(dataCFA$attention_check1)
dataCFA$attentioncheckA <- ifelse(dataCFA$attention_check1 == "right", 0, 1)
table(dataCFA$attentioncheckA)
table(dataCFA$attention_check2)
dataCFA$attentioncheckB <- ifelse(dataCFA$attention_check2 == "right", 0, 1)
table(dataCFA$attentioncheckB)
table(dataCFA$attention_check3)
dataCFA$attentioncheckC<- ifelse(dataCFA$attention_check3 == "right", 0, 1)
table(dataCFA$attentioncheckC)
dataCFA$ac <- dataCFA$attentioncheckA + dataCFA$attentioncheckB + dataCFA$attentioncheckC

library(dplyr)

dataCFA <- subset(dataCFA, ac < 2)
dataCFA <- subset(dataCFA, 
                  trauma_disease == "yes" | trauma_accident == "yes" | trauma_disaster == "yes" |
                  trauma_disaster == "yes" | trauma_explosion == "yes" | trauma_toxic == "yes" |
                  trauma_robbery == "yes" | trauma_death == "yes" | trauma_sexual1 == "yes" |
                  trauma_sexual2 == "yes" | trauma_abuse == "yes" | trauma_child == "yes" |
                  trauma_physical == "yes" | trauma_gun == "yes" | trauma_death2 == "yes" |
                  trauma_witness == "yes" | trauma_details == "yes" |
                  trauma_witness_job == "yes" | trauma_details_job == "yes")
dataCFA <- select(dataCFA, gender, AGE, 
                  PCL1 : PCL20, 
                  QIDS1 : QIDS16, PHQ1 : PHQ9, 
                  GAD1 : GAD7, 
                  risky1 : risky14, 
                  AUDIT1 : AUDIT10)

##Removes rows containting all NAs
dataCFA <- dataCFA[rowSums(is.na(select(dataCFA, PCL1:PCL20))) != ncol(select(dataCFA, PCL1:PCL20)), ]

set.seed(21102020)
split <-  sample(1:nrow(dataCFA), size = round(0.5*nrow(dataCFA)), replace = FALSE)
dataE = dataCFA[split ,]
dataC = dataCFA[-split ,]

write.csv2(dataE, "dataExploratory.csv", row.names = FALSE)
write.csv2(dataC, "dataConfirmatory.csv", row.names = FALSE)

