
############################Data loading############################
## 0. Data loading
cleaned_data <- read.csv("cleaned_data.csv",header = T, encoding = "UTF-8")

scales <- read.csv("scales.csv",header = T, encoding = "UTF-8")

############################1. Behavioural--ACC############################



ACC_Model <- glmer(ACC ~ FT_FB * Target * Indirectness + Implicature + 
                 (1 |Subject)+
                 (1 + FT_FB + Indirectness |Item), family = binomial,
               data =cleaned_data, glmerControl(optimizer = "bobyqa"))
summary(ACC_Model)
confint(ACC_Model, method = "Wald")

mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "0"),"ACC"]))
mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "1"),"ACC"]))

############################2. SCR analyses############################

#### 2.1 SCR latency#### 

summary(lmer(log(SCR_Latency+1) ~ FT_FB * Target * Indirectness + Implicature + 
               (1|Subject),
             data = cleaned_data,control = lmerControl(optimizer = "bobyqa")))

mean(na.omit(cleaned_data$SCR_Latency))

#### 2.2 SCR mean amplitude#### 

EDAMean_model <- lmer(log(SCR_Mean) ~ FT_FB * Target * Indirectness + Implicature + 
                        (1|Subject),
                      data = cleaned_data,control = lmerControl(optimizer = "bobyqa"))

summary(EDAMean_model)
confint(EDAMean_model, method = "Wald")

### 2.2.1 Interaction follow-ups: Target
EDA_sim_Self_FTFB <- lmer(log(SCR_Mean) ~ Indirectness*FT_FB +
                            (1|Subject),
                          data =cleaned_data[which(cleaned_data$Target == "-1"),])
a <- summary(EDA_sim_Self_FTFB)
summary(EDA_sim_Self_FTFB)

EDA_sim_Other_FTFB <- lmer(log(SCR_Mean) ~ Indirectness*FT_FB +
                             (1 |Subject),
                           data =cleaned_data[which(cleaned_data$Target == "1"),])
b <- summary(EDA_sim_Other_FTFB)
summary(EDA_sim_Other_FTFB)

### 2.2.2 Simple effects of Indirectness

EDA_sim_FT_Self <- lmer(log(SCR_Mean) ~ Indirectness + (1|Subject),
                        data =cleaned_data[which(cleaned_data$FT_FB == "-1" &
                                                   cleaned_data$Target == "-1"),])
e <- summary(EDA_sim_FT_Self)
summary(EDA_sim_FT_Self)

EDA_sim_FB_Self <- lmer(log(SCR_Mean) ~ Indirectness + (1|Subject),
                        data =cleaned_data[which(cleaned_data$FT_FB == "1" &
                                                   cleaned_data$Target == "-1"),])
g <- summary(EDA_sim_FB_Self)
summary(EDA_sim_FB_Self)


EDA_sim_FT_Other <- lmer(log(SCR_Mean) ~ Indirectness + (1 |Subject),
                         data =cleaned_data[which(cleaned_data$FT_FB == "-1" &
                                                    cleaned_data$Target == "1"),])
f <- summary(EDA_sim_FT_Other)
summary(EDA_sim_FT_Other)
confint(EDA_sim_FT_Other, method = "profile")


EDA_sim_FB_Other <- lmer(log(SCR_Mean) ~ Indirectness + (1|Subject),
                         data =cleaned_data[which(cleaned_data$FT_FB == "1" &
                                                    cleaned_data$Target == "1"),])
h <- summary(EDA_sim_FB_Other)
summary(EDA_sim_FB_Other)
p.adjust(c(e[[10]][2,"Pr(>|t|)"],
           g[[10]][2,"Pr(>|t|)"],
           f[[10]][2,"Pr(>|t|)"],
           h[[10]][2,"Pr(>|t|)"]), method = "fdr")

### 2.2.3 Descriptive statistics

DesEDA <- describeBy(SCR_Mean  ~ FT_FB * Target * Indirectness, data = cleaned_data,mat=TRUE,digits=4)
colnames(DesEDA) <- c("Condition","FT_FB","Target","Indirectness","n","vars","mean","sd","median","trimmed","mad","min","max","range","skew","kurtosis","se")


### 2.2.4 correlation
SCR_Mean_sub <- aggregate(SCR_Mean ~ Subject * Indirectness * Target * FT_FB, cleaned_data, FUN = mean)
SCR_Mean_sub1<- dcast(SCR_Mean_sub, Subject * FT_FB * Target ~ Indirectness, value.var = "SCR_Mean")
SCR_Mean_sub1$SCR_Mean_diff <- SCR_Mean_sub1$`1` - SCR_Mean_sub1$`0`

SCR_cor <- left_join(SCR_Mean_sub1[,c(1:3,6)],scales, by = "Subject")

SCR_cor_stats <- rbind(cbind(cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==1),'SCR_Mean_diff'],
                             SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==1),'IRI'], method = "pearson")[[3]],
                    cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==1),'SCR_Mean_diff'],
                             SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==1),'IRI'], method = "pearson")[[4]]),
                    cbind(cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                   SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'IRI'], method = "pearson")[[3]],
                          cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                   SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'IRI'], method = "pearson")[[4]]),
                    cbind(cor.test(SCR_cor[which(SCR_cor$FT_FB==1&SCR_cor$Target==1),'SCR_Mean_diff'],
                                   SCR_cor[which(SCR_cor$FT_FB==1&SCR_cor$Target==1),'IRI'], method = "pearson")[[3]],
                          cor.test(SCR_cor[which(SCR_cor$FT_FB==1&SCR_cor$Target==1),'SCR_Mean_diff'],
                                   SCR_cor[which(SCR_cor$FT_FB==1&SCR_cor$Target==1),'IRI'], method = "pearson")[[4]]),
                    cbind(cor.test(SCR_cor[which(SCR_cor$FT_FB==1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                   SCR_cor[which(SCR_cor$FT_FB==1&SCR_cor$Target==-1),'CSF'], method = "pearson")[[3]],
                          cor.test(SCR_cor[which(SCR_cor$FT_FB==1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                   SCR_cor[which(SCR_cor$FT_FB==1&SCR_cor$Target==-1),'CSF'], method = "pearson")[[4]]))

SCR_cor_stats

p.adjust(SCR_cor_stats[,1], method = "fdr")

ols <- lm(SCR_Mean_diff~IRI,data = SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),])
summary(ols)

rrhuber <- rlm(SCR_Mean_diff~IRI,data = SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),])
summary(rrhuber)

sfsmisc::f.robftest(rrhuber, var = "IRI")

SCR_cor_stats_subscale <- rbind(cbind(cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                         SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'PT'], method = "pearson")[[3]],
                                cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                         SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'PT'], method = "pearson")[[4]]),
                                cbind(cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                               SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'FS'], method = "pearson")[[3]],
                                      cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                               SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'FS'], method = "pearson")[[4]]),
                                cbind(cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                               SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'EC'], method = "pearson")[[3]],
                                      cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                               SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'EC'], method = "pearson")[[4]]),
                                cbind(cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                               SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'PD'], method = "pearson")[[3]],
                                      cor.test(SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'SCR_Mean_diff'],
                                               SCR_cor[which(SCR_cor$FT_FB==-1&SCR_cor$Target==-1),'PD'], method = "pearson")[[4]])
                                )
SCR_cor_stats_subscale

############################3. IA1############################


summary(glmer(IA1_SKIP ~ FT_FB * Target * Indirectness + Implicature + 
                (1|Subject)+
                (1|Item), family = binomial,
              data =cleaned_data, glmerControl(optimizer = "bobyqa")))

####3.1 First Fixation####
summary(lmer(log(IA1_FIRST_FIXATION_DURATION) ~ FT_FB * Target * Indirectness + Implicature +
               (1|Subject)+
               (1+ Target|Item),
             data =cleaned_data, control = lmerControl(optimizer = "bobyqa")))

####3.2 gaze####
frlmer <- lmer(log(IA1_GAZE_DURATION) ~ FT_FB * Target * Indirectness + Implicature +
                 (1 |Subject)+
                 (1+Indirectness+Target+ FT_FB|Item),
               data =cleaned_data, control = lmerControl(optimizer = "bobyqa"))

summary(frlmer)
confint(frlmer, method = "Wald")

round(mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "0"),"IA1_GAZE_DURATION"])))

round(mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "1"),"IA1_GAZE_DURATION"])))


####3.3 Second Pass####
srlmer <- lmer(log(IA1_SECOND_PASS_READING_TIME) ~ FT_FB * Target * Indirectness + Implicature +
                 (1+Target |Subject)+
                 (1 +Target |Item),
               data =cleaned_data, control = lmerControl(optimizer = "bobyqa"))

summary(srlmer)

confint(srlmer, method = "Wald")

round(mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "0"),"IA1_SECOND_PASS_READING_TIME"])))

round(mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "1"),"IA1_SECOND_PASS_READING_TIME"])))

####3.4Regression in####

IA1_RI_Model <- glmer(IA1_REGRESSION_IN ~ FT_FB * Target * Indirectness + Implicature +
                        (1 |Subject)+
                        (1 + Target|Item),family = binomial,
                      data =cleaned_data, glmerControl(optimizer = "bobyqa"))

summary(IA1_RI_Model)
confint(IA1_RI_Model, method = "Wald")

##3.4.1 interaction
IA1_RI_sim_Self <- glmer(IA1_REGRESSION_IN ~ Indirectness + Implicature  + (1|Subject)+(1 |Item),
                         data =cleaned_data[which(cleaned_data$Target == "-1"),], family = binomial)

a <- summary(IA1_RI_sim_Self)
summary(IA1_RI_sim_Self)

IA1_RI_sim_Other <- glmer(IA1_REGRESSION_IN ~ Indirectness + Implicature + (1+Indirectness|Subject)+(1 +Indirectness|Item),
                          data =cleaned_data[which(cleaned_data$Target == "1"),], family = binomial)
b <- summary(IA1_RI_sim_Other)
summary(IA1_RI_sim_Other)

confint(IA1_RI_sim_Other, method = "Wald")

p.adjust(c(a[[10]][2,"Pr(>|z|)"],b[[10]][2,"Pr(>|z|)"]), method = "fdr")

mean(na.omit(cleaned_data[which(cleaned_data$Target==1&cleaned_data$Indirectness == "0"),"IA1_REGRESSION_IN"]))
mean(na.omit(cleaned_data[which(cleaned_data$Target==1&cleaned_data$Indirectness == "1"),"IA1_REGRESSION_IN"]))

##3.4.2 correlation

IA1_RI_per_Sub <- aggregate(IA1_REGRESSION_IN ~ Subject*FT_FB*Target*Indirectness, data = cleaned_data, FUN = "mean")
IA1_RI_per_Sub <- aggregate(IA1_REGRESSION_IN ~ Subject*Target*Indirectness, data = IA1_RI_per_Sub, FUN = "mean")

IA1_RI_per_Sub1<- dcast(IA1_RI_per_Sub, Subject * Target ~ Indirectness, value.var = "IA1_REGRESSION_IN")
IA1_RI_per_Sub1$IA1_RI_diff <- IA1_RI_per_Sub1$`1` - IA1_RI_per_Sub1$`0`

IA1_RI_cor <- left_join(IA1_RI_per_Sub1[,c(1,2,5)],scales, by = "Subject")

IA1_RI_cor_stats <- rbind(cbind(cor.test(IA1_RI_cor[which(IA1_RI_cor$Target==1),'IA1_RI_diff'],
                                         IA1_RI_cor[which(IA1_RI_cor$Target==1),'CSF'], method = "pearson")[[3]],
                                cor.test(IA1_RI_cor[which(IA1_RI_cor$Target==1),'IA1_RI_diff'],
                                         IA1_RI_cor[which(IA1_RI_cor$Target==1),'CSF'], method = "pearson")[[4]]),
                       cbind(cor.test(IA1_RI_cor[which(IA1_RI_cor$Target==-1),'IA1_RI_diff'],
                                      IA1_RI_cor[which(IA1_RI_cor$Target==-1),'CSF'], method = "pearson")[[3]],
                             cor.test(IA1_RI_cor[which(IA1_RI_cor$Target==-1),'IA1_RI_diff'],
                                      IA1_RI_cor[which(IA1_RI_cor$Target==-1),'CSF'], method = "pearson")[[4]]))

IA1_RI_cor_stats

p.adjust(SCR_cor_stats[,1], method = "fdr")

####Total Reading####
IA1_TR_Model <- lmer(log(IA1_TOTAL_READING_TIME) ~ FT_FB * Target * Indirectness + Implicature +
                       (1 + Indirectness|Subject)+
                       (1 + Indirectness|Item),
                     data =cleaned_data, control = lmerControl(optimizer = "bobyqa"))
summary(IA1_TR_Model)
confint(IA1_TR_Model, method = "Wald")

round(mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "0"),"IA1_TOTAL_READING_TIME"])))

round(mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "1"),"IA1_TOTAL_READING_TIME"])))


############################ 4. IA2############################


summary(glmer(IA2_SKIP ~ FT_FB * Target * Indirectness + Implicature +
                (1|Subject)+
                (1|Item), family = binomial,
              data =cleaned_data, glmerControl(optimizer = "bobyqa")))

####4.1 First Fixation####
summary(lmer(log(IA2_FIRST_FIXATION_DURATION) ~ FT_FB * Target * Indirectness + Implicature +
               (1+FT_FB|Subject)+
               (1|Item),
             data =cleaned_data, control = lmerControl(optimizer = "bobyqa")))


####4.2 Gaze####
IA2_GD_Model <- lmer(log(IA2_GAZE_DURATION) ~ FT_FB * Target * Indirectness + Implicature +
                       (1 + Indirectness|Subject)+
                       (1 + FT_FB + Indirectness |Item),
                     data =cleaned_data, control = lmerControl(optimizer = "bobyqa"))
summary(IA2_GD_Model)

confint(IA2_GD_Model, method = "Wald")

round(mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "0"),"IA2_GAZE_DURATION"])))

round(mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "1"),"IA2_GAZE_DURATION"])))

####4.3 Second Pass####
summary(lmer(log(IA2_SECOND_PASS_READING_TIME) ~ FT_FB * Target * Indirectness + Implicature +
               (1 |Subject)+
               (1 |Item),
             data =cleaned_data, control = lmerControl(optimizer = "bobyqa")))


####4.4 Regression in####
summary(glmer(IA2_REGRESSION_IN ~ FT_FB * Target * Indirectness + Implicature + 
                (1+Target|Subject)+
                (1+FT_FB|Item), family = binomial,
              data =cleaned_data, glmerControl(optimizer = "bobyqa")))

####4.5 Regression out####
IA2_RO_Model <- glmer(IA2_REGRESSION_OUT ~ FT_FB * Target * Indirectness + Implicature + 
                        (1|Subject)+
                        (1+Indirectness|Item), family = binomial,
                      data =cleaned_data, glmerControl(optimizer = "bobyqa"))
summary(IA2_RO_Model)
confint(IA2_RO_Model, method = "Wald")

mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "0"),"IA2_REGRESSION_OUT"]))
mean(na.omit(cleaned_data[which(cleaned_data$Indirectness == "1"),"IA2_REGRESSION_OUT"]))


####4.6 Total Reading####
summary(lmer(log(IA2_TOTAL_READING_TIME) ~ FT_FB * Target * Indirectness + Implicature +
               (1+Indirectness |Subject)+
               (1+Indirectness  |Item),
             data =cleaned_data, control = lmerControl(optimizer = "bobyqa")))
