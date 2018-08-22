##### PACKAGES #####
# Put all the necessary packages here
# So that other could easily install them before compiling multiple times
# I prefer to use require() insted of library() because
# it will only load the library if it hasn't already.
require(ggplot2)
require(magrittr)

##### DATA LOAD ######
# Load your data in this section - it is crucial to know on what data you are actually working
#load("analysis/data/raw_data/data3.R")
# New database
load("analysis/data/raw_data/data4.R")
data3 <- data4

##### FUNCTIONS #####
# Put all the functions above the code below,
# so that you can easily find them later and execute the code seemlesly

sum_adrenalines <- function(vars,data){
  lapply(data[,vars],function(var){
    ifelse(!is.na(var),
           ifelse(var=="yes", 1, 0),
           0)}) %>%
    Reduce(f=`+`) %>%
    return()
}

rcalc <- function(var,level="yes"){
  tab <- table(data3[,var],data3$rANA)
  if(level %in% (data3[,var] %>% levels)){
    pn <- tab[level,"no"]
    pp <- tab[level,"yes"]
    l <- level
  } else {
    pn <- tab[2,"no"]
    pp <- tab[2,"yes"]
    l <- levels(data3[,var])[2]
  }
  nn <- data3$rANA %>% factor %>% summary %>% {.[1]}-pn
  np <- data3$rANA %>% factor %>% summary %>% {.[2]}-pp
  eval <- c(pn,pp,nn,np) %>% matrix(byrow = T,ncol = 2)
  prt <- (prop.table(eval,2)*100) %>% {signif(.,3)}
  p <- eval %>% fisher.test() %>% .$p.value %>% signif(3)
  return(c(var,prt[1,1],prt[1,2],p,l))
}

f1 <- function(x){
  names(x) <- NULL
  return(do.call(t.test,x)$p.value %>% signif(3))
}
f2 <- function(x){
  x %>% table(data3$rANA) %>% {.[1:2,]} %>% fisher.test() %>% {.$p.value} %>% signif(3)
}
f3 <- function(x){
  x %>%
    lapply(function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
}


##### VARIABLES #####
# If you are making new variables than put them here for reference later.
# It is always good to have all derived data show up early in the script
# so that you would not have issues with code execution later.

adr2 <- which(data3$d_560_adren2_v5=="yes") # Which patients recieved 2 times the adrenilin?
data.r <- data3[adr2,]
#table(data.r$d_severity_rm)

sev <- data.r %>% {which(.$d_severity_rm%in%c("grade III","grade IV"))}
data.r <- data.r[sev,]


# So First get the idea who got 2 shots of adrenaline or more
adre.no <- sum_adrenalines(vars= c("q_521_autoinj_v5",    # Patient self medicated
                                   "q_522_adren_iv",   # First line treatent by prof
                                   "q_522_adren_im",      #
                                   "q_522_adren_inhal",
                                   "q_530_adren2_v5",     # seond dose of adrenalin ?
                                   "q_552_adren_im_v5",   # second line treatemt
                                   "q_552_adren_iv_v5",
                                   "q_552_adren_inhal_v5"
                        ),
                data = data3)

data3$d_adre_no <- adre.no
sumAdre <- summary(data3$d_adre_no %>% factor())

# Then figure out how close were these to one another in time
# data2adr <- data3[reps <- c(which(!is.na(data3$q_530_adren2_time_in_min_v5)),
                            # which(adre.no>1),
                            # which(data3$q_530_adren2_v5=="yes")) %>% unique(),]
#late <- which(data$q_530_adren2_time_in_min_v5>20) # exclude these late adrinalin shots
#data2adr[late,"q_140_fatal"]
#data2adr[late[1],]
#data2adr[late,
#         refractory <- data2adr$d_adre_no>2| # sure! diagnosis of refractory ana
#              data2adr$q_140_fatal=="yes"|
#              ifelse(!is.na(data2adr$q_530_adren2_time_in_min_v5)
 #                    data2adr$q_530_adren2_time_in_min_v5<20,
  #                   FALSE)
# next see if the patients had a severe response III or IV
# This will be our cohort
# Next determine the reaction based on case estimation by evaluating each field visually.
# next see the fatalities - this group is going to be of interest !

#data3$d_520_adren1 %>% summary
#data3$d_522_adren_agg %>% summary

#Example:
varSym <- c(365:367,17:70)
varTher <- c(227:299)

# add Seveirity
# Define a new variable - reaction type to exclude non ANA cases
data <- data3
sentinel_skin <- rep(F, nrow(data))
sentinel_skin[which(data$q_111_angioedema=="yes"|
                      data$q_111_urticaria=="yes"|
                      data$q_111_erythema_flush_v5=="yes"|
                      data$q_111_pruritus=="yes")] <- T

sentinel_respiratory <- rep(F, nrow(data))
sentinel_respiratory[which(data$q_115_cyanosis_pallor_v6=="yes"|
                             data$q_113_stridor_inspiratory=="yes"|
                             data$q_113_wheezing_expiratory_distress_v5=="yes"|
                             data$q_113_dyspnea=="yes"|
                             data$q_113_respiratory_arrest=="yes")] <- T

sentinel_cardio <- rep(F, nrow(data))
sentinel_cardio[which(
  data$q_112_incontinence=="yes"|
    data$q_114_reductions_of_alertness=="yes"|
    data$q_114_loss_of_consciousness=="yes"|
    data$q_114_cardiac_arrest=="yes"|
    data$q_114_hypotension_collapse_v5=="yes")] <- T

sentinel_gastro <- rep(F, nrow(data))
sentinel_gastro[which(data$q_112_abdominal_pain=="yes"|
                        data$q_112_vomiting=="yes")] <- T

reaction_type_brown <- rep(NA, length(data[,1]))
# Correct the sum of organs involved.
organ_sum <- rep(0,nrow(data))
organ_sum <- organ_sum+ ifelse(!is.na(data$q_111),
                               ifelse(data$q_111=="yes",1,0),
                               0)
organ_sum <- organ_sum+ ifelse(!is.na(data$q_112),
                               ifelse(data$q_112=="yes",1,0),
                               0)
organ_sum <- organ_sum+ ifelse(!is.na(data$q_113),
                               ifelse(data$q_113=="yes",1,0),
                               0)
organ_sum <- organ_sum+ ifelse(!is.na(data$q_114),
                               ifelse(data$q_114=="yes",1,0),
                               0)

organ_sum %<>% factor()
reaction_type_brown[which(organ_sum=="1"&data$d_111_sum%in%
                            c("one","two","three","four","five"))] <- "skin only"

reaction_type_brown[which(sentinel_skin==T&(sentinel_cardio==T|sentinel_respiratory==T))] <- "anaphylaxis"

data$q_120_time_between_v4 <- factor(data$q_120_time_between_v4,
                                     levels = c("-9","0","1","2","3","4","5","6","7"),
                                     labels=c("NA","unknown","10","30",
                                              "60","120","4h","more","4h-5"))
levels(data$q_120_time_between_v4) <- levels(data$q_120_time_between_v4)[c(1:8,7)]

# BELOW: I am not including the information about the time from exposure to reaction as it was largly missing. The cases are nevertheless supposed allergic reactions therefore the likelyhood of an allergen and occurence of these symptoms should be strict enough to be sufficient as a definition of anaphylaxis.

reaction_type_brown[which(#data$q_120_time_between_v4%in%levels(data$q_120_time_between_v4)[3:7]&
  ((sentinel_skin==T&sentinel_cardio==T)|
     (sentinel_skin==T&sentinel_gastro==T)|
     (sentinel_skin==T&sentinel_respiratory==T)|
     (sentinel_cardio==T&sentinel_gastro==T)|
     (sentinel_cardio==T&sentinel_respiratory==T)|
     (sentinel_respiratory==T&sentinel_gastro==T)|
     (data$q_140_fatal=="yes")))] <- "anaphylaxis"
reaction_type_brown[is.na(reaction_type_brown)] <-"ANA-def not met"

reaction_type_brown <- as.factor(reaction_type_brown)
data$reaction_type_brown <- reaction_type_brown

# Add severity of anaphylaxis
severity_brown <- rep(NA,length(data[,1]))
severity_brown[which(reaction_type_brown=="anaphylaxis")]<-"mild"
severity_brown[which(reaction_type_brown=="anaphylaxis"&
                       (data$q_115_cyanosis_pallor_v6=="yes"|
                          data$q_114_hypotension_collapse_v5=="yes"|
                          data$q_114_loss_of_consciousness=="yes"|
                          data$q_114_reductions_of_alertness=="yes"|
                          data$q_112_incontinence=="yes"|
                          data$q_140_fatal=="yes"|
                          data$q_114_cardiac_arrest=="yes"|
                          data$q_113_respiratory_arrest=="yes"))] <- "severe"
severity_brown <- as.factor(severity_brown)
data$severity_brown <- severity_brown
data3 <- data




#data3[6515,varSym]
#data3[6515,varTher] # multiple doses of adrenalin and other treatment. she recovered. refractory!
#data3[2925,varSym]
#data3[2925,varTher] # - food case fatal !  although two doses administered by patient and professional.

#data3[3582,varSym]
#data3[3582,varTher] # - food case fatal !  although two doses administered by patient and professional.


# Patient has to get at least 2 doses of adrenalin and the reaction has to be either fatal or
# refractory anaphylaxis is when you give adrenaline and the patient dies
# Helpful variables in our case
#"q_521_autoinj_v5" %>% rcalc
#"q_521_other_v5" %>% rcalc
data3$q_521_other_free_v5 %>% summary # Mund-zu-Mund-Beatmung, ABC
data3$q_522_adren_im %>% summary
data3$q_522_adren_iv %>% summary
data3$q_522_adren_inhal %>% summary
data3$q_522_volume %>% summary
data3$q_522_other_free %>% summary # insgesamt 3 Adrenalin-injektionen durch Hausarzt,
# Ephedrine ; Ephedrin, Dosis unbekannt; Ephedrin 9mg i.v., Infusionen
# Ephedrin 50mg, Infusionen,Adrenalin als wiederholter Bolus (100 mikrogramm)ok, dann Adrenalinperfusor, Noradrenalinperfusor, Trendelenburglagerung
# Ephedrin 40mg, Magnesium 2 Ampullen, Defibrillation 3x 200Joule, Amiodaron 1 Ampulle
# Ephedrin
# Dobutamin im Rahmen des Koronarsyndroms
# Atropin mehrfach 0,5 mg
# Atropin i.v. , Spasfon i.v. (Phloroglucinol), Metoclopramid i.v
# Atropin Boli, 2 minütige mechanische und medikamentöse Reanimation
# Atropin (Spasmolytikum)
# Adrenalin, Ranitic, Atropin ok
# Adrenalin s.c., Buscopan i.v.
# Adrenalin i.m.
# Adrenalin (verabreichungsmodus unbekannt), Ulsal, Paspertin, Respicur
# Adrenalin
# ephedrine
# Reanimation, weitere Medi unbekannt
# Atropin
# Adrenalin s.c.
data3$q_530_adren2_v5 %>% summary # Very important!!!
data3$q_530_adren2_who_v5 %>% summary
data3$q_530_adren2_time_in_min_v5 %>% summary # This can tell us if the effect was satisfactory !
data3$q_550_2nd_v5 %>% summary # second line treatment
data3$q_552_adren_im_v5 %>% summary
data3$q_552_adren_iv_v5 %>% summary
data3$q_552_adren_inhal_v5 %>% summary
data3$q_561_hospital_admission_v6 %>% summary
data3$q_562_intensive_care_v6 %>% summary

data3$q_comments_v5 %>% summary()
data3$Q_COM0 %>% summary()
data3$Q_COM2 %>% summary
data3$Q_COM1 %>% summary
data3$d_552_adren_agg_v5 %>% summary
data3$q_530_adren2_time_in_min_v5 %>% hist(breaks=1000)
which(data3$q_530_adren2_time_in_min_v5>20) # get rid of these that are over 10 minutes apart as these probably responded well.
data3$q_143_fatal_adre_v5 %>% summary
data3$q_140_fatal %>% summary

export <- list()

refractory.death <- table(data3$d_adre_no,data3$q_140_fatal)
export$refractoryDeathTab <- refractory.death
#which(data3$d_adre_no==2& data3$q_140_fatal=="yes")

# refractory anaphylaxis is when you give adrenalin 2 times and there has been a second line treatment with
#some meds other than adrenalin - iv antihist,
export$drugsRefractory2ndLine <- c("q_552_antih_iv_v5",
                                          "q_552_beta2_iv_v5",
                                          "q_552_dopamine_v5",
                                          "q_552_glucagon_v5",
                                          "q_552_cortico_iv_v5",
                                          "q_552_methyleneb_v5",
                                          "q_552_theo_iv_v5",
                                          "q_552_volume_v5")
# ! Warning non
#czeck which cases had at least one of these drugs in their secon line treatment.
data3$drug2ndLine <- sum_adrenalines(export$drugsRefractory2ndLine,data = data3)
export$Tab2Line <- table(Drugs2Line = data3$drug2ndLine,
      Adrenalin = data3$d_adre_no)

export$Tab2Line[3:4,3:5] %>% sum # Lets look into these patients
which(data3$drug2ndLine==2&data3$d_adre_no==2)
data3[3049,varSym]
data3[3049,varTher]
# This patient got 2 Adrenalins + cortico iv + antihistaminics iv in a second line treatment. not refractory.
data3[3099,varSym]
data3[3099,varTher]
# This patient got 2 Adrenalins im + inhal + cortico iv +volume in a second line treatment. not refractory. questionable GRADE III
data3[4196,varSym]
data3[4196,varTher]
# This patient got 2 Adrenalins + cortico iv + antihistaminics iv in a second line treatment. not refractory.
which(data3$drug2ndLine==3&data3$d_adre_no==2)
data3[9893,varSym]
data3[9893,varTher] # This might be refractory but... grade III got voluem +beta2, antih. addmitted to the hospital adren im second dose


# Or we state that the patients needed to be addmited to the hospital if did not die?
export$hospitalTab <- table(Hospital=data3$q_561_hospital_admission_v6,
                            Adrenalin_Doses = data3$d_adre_no)

which(data3$q_561_hospital_admission_v6=="yes"&data3$d_adre_no>=2)
data3[2693,varSym]
data3[2693,varTher] # not a good example as the second line treatment is lacking here

data3[10108,varSym]
data3[10108,varTher] # Grade II reaction ...


# We either say we did second line + hospital admission or we say we had reactions type IV
table(Severity = data3$d_severity_rm, data3$d_adre_no)
which(data3$d_severity_rm=="grade IV"& data3$d_adre_no==4) # This patient was already evaluated above - this is not refractory.
which(data3$d_severity_rm=="grade IV"& data3$d_adre_no==3)
data3[9771,varTher] # Adrenalin only. non- refractory

export$gradeIV2drugs <- table(Drug2ndLine = data3$drug2ndLine[data3$d_severity_rm=="grade IV"],
      AdreDoses = data3$d_adre_no[data3$d_severity_rm=="grade IV"])

refractoryCasesRM <- c(which(data3$d_adre_no>1& data3$q_140_fatal=="yes"),
                     which(data3$d_adre_no>1&data3$drug2ndLine>1&(data3$d_severity_rm=="grade IV"|
                                                 (data3$d_severity_rm=="grade III"&data3$q_561_hospital_admission_v6=="yes"))))


data3[5694,varSym]
data3[5694,varTher] # This guy got only inhalaitve and im adrenalin and grade III - questionable.

data3[5884,varSym]
data3[5884,varTher] # This guy got 3 doses of adrenaline and grade III and emergency treatemnt.


varsImp <- c(11,12,366,367,368,42,43,47,55,56,57,261,262,264,271,297,298,299,66)

refractoryCasesBrown <- c(which(data3$d_adre_no>1& data3$q_140_fatal=="yes"),
                          which(data3$d_adre_no>1&data3$drug2ndLine>1&data3$severity_brown=="severe"))
export$cases <-  data3[refractoryCasesBrown,varsImp]

data3[2696,varSym]
data3[2696,varTher] # This is refractory. multiple adrenalin no info about hospital treatemnt.

data3[3049,varSym]
data3[3049,varTher] # This is not refractory. multiple adrenalin no info about hospital treatemnt, glucocorticoids and antihistamines.

refractoryCasesBrown%in%refractoryCasesRM
refractoryCasesRM%in%refractoryCasesBrown

export$casesBrown <- data3[refractoryCasesBrown,]
save(export, file = "analysis/data/derived_data/export.R")




##### Perioperative variable ####
perioperative <- rep("no",length(data3[,1]))
perioperative[which(data3$q_152_location=="medical practice, hospital"&
                (!is.na(data3$q_332_analgesic)|
                   data3$q_3362_gen_anaesthetics_v5=="yes"|
                      data3$q_3361_local_anaesthetics_v4=="yes"))] <- "yes"
perioperative %>% factor() %>% summary()
data3$perioperative <- perioperative %>% factor()

# Other tests.
# rdb$q_210_diagnosis_v5 %>% summary
# data3$q_552_methyleneb_v5 %>% summary()
# rdb$q_522_glucagon_v5 %>% summary()
# rdb[which(rdb$q_530_adren2_time_in_min_v5 >15),varTher]
# rdb$q_152_location %>% summary() #=="medical practice, hospital"
#
#
#
# rdb$q_116_VAS_v7 %>% summary()
#
# rdb$d_620_autoinj_total
# rdb$q_521_autoinj_v5
# rdb$q_620_autoinj_acute
# rdb$q_620_autoinj_prior_v5
# rdb$q_632_autoinj_number_v7
# rdb$q_2311_who_diagnosed_v5 %>% summary()
# rdb$d_adre_no %>% factor %>% summary()
# rdb$q_212_tryptase_value_v5 %>% plot(x = rdb$d_severity_rm)
# rdb$d_210_tryptase_cat %>% summary()
#
#

########### Definition of Refractory anaphylaxis ##############
# 1. Death and two adrenalin doses.
refractory.death <- rep("no",length(data3[,1]))
refractory.death[data3$d_adre_no>=2 & data3$q_140_fatal == "yes"] <- "yes"

# 2. 2 or more Epi, severe reaction, hospital admission
refractory <-rep("no",length(data3[,1]))
refractory[data3$d_adre_no>=2 &
             data3$severity_brown =="severe" &
             data3$q_561_hospital_admission_v6 == "yes"] <- "yes"
summary(refractory %>% factor)
table(refractory,refractory.death)

# 3. severe + refractory 2ndline medication
refractory.med <- rep("no",length(data3[,1]))
refractory.med[data3$d_adre_no>1&data3$drug2ndLine>1&data3$severity_brown=="severe"] <- "yes"
table(refractory,refractory.med)

# 4. 2 or more Epi, severe reaction, hospital admission + 2nd line Drugs
refractory.strict <-rep("no",length(data3[,1]))
refractory.strict[data3$d_adre_no>=2 &
             data3$severity_brown =="severe" &
             data3$q_561_hospital_admission_v6 == "yes"&
             data3$drug2ndLine>1] <- "yes"
summary(refractory.strict %>% factor)
table(refractory.strict,refractory.death)

# 5. 2 or more Epi, severe reaction, hospital admission + Intensive care
refractory.ICU <-rep("no",length(data3[,1]))
refractory.ICU[data3$d_adre_no>=2 &
                    data3$severity_brown =="severe" &
                    data3$q_562_intensive_care_v6 == "yes"] <- "yes"
summary(refractory.ICU %>% factor)
table(refractory.strict,refractory.ICU)

table(refractory.ICU,refractory.death)

# Sum up the fatalities and strict refractory cases
refractory.f <- rep("no",length(data3[,1]))
refractory.f[refractory.strict=="yes"|refractory.death=="yes"|refractory.ICU=="yes"] <- "yes"
summary(refractory.f %>% factor)

data3$q_comments_v5[refractory.f=="yes"]
data3$Q_COM0 [refractory.f=="yes"]
data3$Q_COM2 [refractory.f=="yes"]
data3$Q_COM1 [refractory.f=="yes"]


###### Manual selection process #################
rdb <- data3[refractory.f=="yes",]
#rdb[rdb$b_patient_code=="m-1963-g-10-e-08-31 ",] # The Adrenaline dosis could not be given because the needle broke...
#rdb[rdb$b_patient_code=="25c-2012 (88c)      ",] # Emergency treatment solely professional adrenalin was given later.
#rdb[rdb$b_patient_code=="24c-2012 (87c)      ",] # Only emergency treatment also ... no autoinjector

#rdb[rdb$d_520_adren1=="no",][1,] #<- This is rather not refractory as the patient recieved adrenalin and inhalatory adrenalin...
#rdb[rdb$d_520_adren1=="no",][2,] # Non refractory ! Multiple reactions to potentially two allergens.
#rdb[rdb$d_520_adren1=="no",][3,] # Refractory! 3 varoius routes of administration... Peanuts. Delayed Adrenalin administration

#rdb[rdb$d_adre_no==2,][c(7,8),]

#rdb[rdb$q_130_biphasic_v4=="yes",][c(1,2),] # -1 ADREN INHALATIV?!?!?!
#rdb[rdb$q_130_biphasic_v4=="yes",][c(5,6),] #
# 1. pat - multiple triggers,
# biphasic reactions that lead to multiple adrenalin doses but were effective
#
REMOVE <- rdb$b_case_id[which(rdb$q_130_biphasic_v4=="yes"&rdb$q_140_fatal=="no")]
R1 <- REMOVE %>% length()
rdb <- rdb[!rdb$b_case_id%in%REMOVE,]
#rdb[which(rdb$q_530_adren2_time_in_min_v5>=30),c("b_case_id","d_560_adren2_v5","q_551_2nd_who_v5","q_152_location","d_severity_rm","d_elicitor_gr5","q_530_adren2_v5","q_130_biphasic_v4","d_522_adren_agg","q_552_volume_v5","d_552_adren_agg_v5","q_562_intensive_care_v6","drug2ndLine","q_522_adren_im","q_522_adren_iv","q_550_2nd_v5","q_140_fatal","q_530_adren2_time_in_min_v5","d_adre_no")] # Both non refractory
#rdb[which(rdb$q_530_adren2_time_in_min_v5>=30),][1,]
#rdb[rdb$q_530_adren2_time_in_min_v5>=30,][c(8,15),] # Both non refractory
#rdb[rdb$q_530_adren2_time_in_min_v5>=30,][c(16,17),] # 15 non ref, 17 REF!
#rdb[rdb$q_530_adren2_time_in_min_v5>=30,][c(18,19),] #
#rdb[rdb$b_case_id=="9841",]

#REMOVE <- c(7905,11538,12813,14260,12331) # Long gap between adrenalin doses
REMOVE <- rdb$b_case_id[which(rdb$q_530_adren2_time_in_min_v5>=10&rdb$q_140_fatal=="no")]
R2 <- length(REMOVE)
rdb <- rdb[!rdb$b_case_id%in%REMOVE,]
REMOVE <- rdb$b_case_id[which(!rdb$q_140_fatal=="yes"& ### Patient responded adequately after two doses of Adreanline
                                rdb$d_560_adren2_v5=="no")]
R3 <- REMOVE %>% length
###### FInal database for further evaluation #######
rdb <-rdb[!rdb$b_case_id%in%REMOVE,]
# rdb[!rdb$q_140_fatal=="yes"&
      # rdb$q_550_2nd_v5=="no",c("b_case_id","d_elicitor_gr5","q_152_location","d_severity_rm","q_562_intensive_care_v6","drug2ndLine",
                              # "d_522_adren_agg","q_522_adren_im","q_522_adren_iv","q_530_adren2_v5","q_550_2nd_v5","d_552_adren_agg_v5","d_560_adren2_v5",
                              # "q_551_2nd_who_v5","q_552_volume_v5","q_140_fatal","q_530_adren2_time_in_min_v5","d_adre_no")]


data3$rANA <- rep (NA,length(data3$b_submitdate))
data3$rANA[data3$severity_brown=="severe"] <- "no"
data3$rANA[data3$b_case_id %in% rdb$b_case_id] <- "yes"

# And the control group
control <- data3[!(data3$b_case_id %in% rdb$b_case_id),]



##### EXPLORATION #####
# Here you may start with a really messy analysis that will be only internally seen.
# Get the feel of the data and identify interesting approaches



# ** Results ----
# Put all the values that you want to export to the paper here.
# All necessary variables and values should be put here for reference
# and to allow easier verification of the results

export <- list()
export$adrenalinNumber <- adre.no %>% as.character %>% factor %>% summary()
export$severityTab <- table(data3$d_adre_no,data3$d_severity_rm)

v<- list()
v$allCases <- length(data3$b_submitdate)
v$casesR <- length(which(data3$rANA=="yes"))
v$casesC <- length(which(data3$rANA=="no"))
v$perioperative <- rcalc("perioperative")
v$casesInsect <- table(data3$q_340_insects[data3$rANA=="yes"])
v$mortality <- rcalc("q_140_fatal")# ** Tables ----
# This section deals with tables. Prepare them here and reference to them in the paper.Rmd

sdb <- data3[which(data3$rANA=="no"),]

#### 1. Demography ################
demoTab <- cbind(n = rdb$b_sex %>% summary(),
                Age = rdb$d_age %>% split(rdb$b_sex) %>%
                              lapply(.,function(x){mean(x) %>% signif(3)}),
                Cardiologic = rdb$q_410_cardio_cur%>% split(rdb$b_sex) %>%
                              f3,
                DM = rdb$q_410_diab_cur_v6%>% split(rdb$b_sex) %>%
                              f3,
                `Food allergy` = rdb$q_410_foodallergy_cur_v6%>% split(rdb$b_sex) %>%
                              f3,
                 Mastocytosis = rdb$q_410_masto_cur %>% split(rdb$b_sex) %>% f3,
                 Malignancy = rdb$q_410_malig_cur %>% split(rdb$b_sex)  %>% f3,
                `Atopic dermatitis` = rdb$q_410_ad_cur %>% split(rdb$b_sex)  %>% f3
                )


demoTabs <- cbind(n = sdb$b_sex %>% summary(),
                 Age = sdb$d_age %>% split(sdb$b_sex) %>%
                   lapply(.,function(x){mean(x) %>% signif(3)}),
                 Cardiologic = sdb$q_410_cardio_cur%>%  split(sdb$b_sex) %>%f3,
                 DM = sdb$q_410_diab_cur_v6%>% split(sdb$b_sex) %>% f3,
                 `Food allergy` = sdb$q_410_foodallergy_cur_v6 %>%  split(sdb$b_sex) %>%f3,
                 Mastocytosis = sdb$q_410_masto_cur %>% split(sdb$b_sex) %>% f3,
                 Malignancy = sdb$q_410_malig_cur %>% split(sdb$b_sex)  %>% f3,
                 `Atopic dermatitis` = sdb$q_410_ad_cur %>% split(sdb$b_sex)  %>% f3
)

demoTabsP <- cbind(n = data3$b_sex %>% table(data3$rANA) %>% summary() %>% {.$p.value} %>% signif(3),
                  Age = data3$d_age %>% split(data3$rANA) %>% f1,
                  Cardiologic = data3$q_410_cardio_cur %>% f2,
                  DM = data3$q_410_diab_cur_v6 %>% f2,
                  `Food allergy` = data3$q_410_foodallergy_cur_v6 %>% f2,
                  Mastocytosis = data3$q_410_masto_cur %>% f2,
                  Malignancy = data3$q_410_malig_cur %>% f2,
                  `Atopic dermatitis` = data3$q_410_ad_cur %>% f2

)
demoTab <- rbind(demoTab,demoTabs,demoTabsP)
demoTab <- cbind(Group = c("refractory","refractory","severe","severe","p value"),demoTab)

rdb$q_410_masto_cur%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
rdb$q_410_asthma_cur%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
rdb$q_410_rhinitis_cur%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
rdb$q_410_ad_cur%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
rdb$q_410_foodallergy_cur_v6%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
rdb$q_410_urtic_cur%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
rdb$q_410_infect_cur%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
rdb$q_410_thyroid_cur%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
rdb$q_410_malig_cur%>% split(rdb$b_sex) %>% lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})

countries <- rdb$d_centres_country %>% summary() %>% .[{which(.!=0)}]

#rdb[rdb$d_severity_rm=="grade II",]



#### 2. Elicitor tab ####

elicitorTab <- cbind(n = rdb$d_elicitor_gr5 %>% summary(),
                     percent = rdb$d_elicitor_gr5 %>% {summary(.)/42*100} %>% signif(3),
                     percANA = control$d_elicitor_gr5 %>% {summary(.)/length(.)*100} %>% signif(3),
                     p = c(rcalc("d_elicitor_gr5","food")[4],
                           rcalc("d_elicitor_gr5","drugs")[4],
                           rcalc("d_elicitor_gr5","insects")[4],
                           rcalc("d_elicitor_gr5","other")[4],
                           rcalc("d_elicitor_gr5","unknown")[4]),
                 Age = rdb$d_age %>% split(rdb$d_elicitor_gr5) %>%
                   lapply(.,function(x){mean(x) %>% signif(3)}),
                 Male = rdb$b_sex%>% split(rdb$d_elicitor_gr5) %>%
                   lapply(.,function(x){(length(which(x=="male"))/length(x)*100) %>% signif(3)}),
                 Peroperative = rdb$perioperative%>% split(rdb$d_elicitor_gr5) %>%
                   lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)}),
                 `Food allergy` = rdb$q_410_foodallergy_cur_v6%>% split(rdb$d_elicitor_gr5) %>%
                   lapply(.,function(x){(length(which(x=="yes"))/length(x)*100) %>% signif(3)})
)

#### 3. Risk Factors #####

# Patient Dependent
# Coronary artery disease


riskTab <- rbind(rcalc("q_410_cardio_cur"),
rcalc("d_age_gr5b","seniors 65+"), # 65+
rcalc("q_410_asthma_cur"),
rcalc("q_423_beta"),
rcalc("perioperative"),
rcalc("d_elicitor_gr5")
)

#### 4. Therapy tab ####
therapyTab <- rbind(
  rcalc("q_522_adren_im"),
  rcalc("q_522_adren_iv"),
  rcalc("q_552_adren_iv_v5"),
  rcalc("q_522_volume"),
  rcalc("q_552_volume_v5"),
  rcalc("q_522_antih_iv"),
  rcalc("q_552_antih_iv_v5"),
  rcalc("q_521_cortic_v5"),
  rcalc("q_522_cortico_iv"),
  rcalc("q_552_cortico_iv_v5"),
  rcalc("q_522_beta2_iv"),
  rcalc("q_552_beta2_inhal_v5"),
  rcalc("q_522_theo_iv"),
  rcalc("q_522_o2"),
  rcalc("q_552_dopamine_v5"),
  rcalc("q_552_glucagon_v5"),
  rcalc("q_552_methyleneb_v5"),
  rcalc("q_561_hospital_admission_v6"),
  rcalc("q_562_intensive_care_v6")
)

therapyTab[,1] <-
  c("adrenaline im.",
    "adrenaline iv.",
    "adrenaline iv., 2nd line",
    "volume",
    "volume, 2nd line",
    "antihistaminics iv.",
    "antihistaminics iv. 2nd line",
    "corticosteoids, all routes",
    "corticosteroids iv.",
    "corticosteroids iv. 2nd line",
    "beta-2-mimetics iv.",
    "beta-2-mimetics inh. 2nd line",
    "theophyline iv.",
    "100% oxygen",
    "dopamine iv.",
    "glucagon iv.",
    "methylene blue",
    "hospital admission",
    "intensive care")

#### 5. Fatal cases ####
data3$q_140_fatal %>% table(data3$rANA)

#### 6. Exact elicitor ####
elicitExact <- rbind(
  rcalc("d_330_drug_group","antibiotics"),
  rcalc("d_330_drug_group","xray_cm"),
  rcalc("d_330_drug_group","muscle relaxant"),
  rcalc("d_320_food_group", "legumes"),
  rcalc("q_340_insects", "bee"),
  rcalc("q_340_insects", "yellow jacket")
)

#### 7. Symptoms tab ####
sympt.all <- names(data3)[18:70] %>% lapply(rcalc) %>% do.call(what = rbind)
symptTab <- sympt.all[c(4,6,14,16,18,21,22,24,26,29,30,31),] %>% rbind(rcalc("q_140_fatal"))
symptTab[,1] <- c("Pruritus",
                  "Skin symptoms",
                  "Respiratory symptoms",
                  "Respiratory arrest",
                  "Chest tightness",
                  "Throat tightness",
                  "Expiratory distress",
                  "Inspiratory stridor",
                  "Loss of consciusness",
                  "Cardiac arrythmia",
                  "Cardiac arrest",
                  "Vertigo",
                  "Death")

# 8. Cofactors tab ####
cofactorsTab <-  names(data3)[c(185:213,215:232)] %>% lapply(rcalc) %>% do.call(what = rbind)
cofactorsTabF <- cofactorsTab[c(-1,-2,-3,-4,-5,-7,-9,-11,-13,-16,-17,-19,-20,-21,-22,-23,-24,-25,-26,-27,
                                -29,-31,-34,-36,-38,-39,-41,-42,-44,-45,-46),]
cofactorsTabF[,1] <- c(
                       "Concomitant asthma",
                       "Concomitant AD",
                       "Concomitant diabetes",
                       "Concomitant cardiologic condition",
                       "Concomitant infection",
                       "History of malignant disease",
                       "Concomitant mastocytosis",
                       "Concomitant other disease - unspecified",
                       "Exercise prior to reaction",
                       "Psychological burden",
                       "Concomitant medication",
                       "ASA",
                       "Beta-blockers",
                       "PPI",
                       "Other drugs",
                       "Alcohol use prior to the reaction")
# ** Figures ----
library(DiagrammeRsvg)
require(DiagrammeR)



# Create a node data frame
nodes <-
  create_node_df(
    n = 8,
    type = "a",
    label= c(paste0("All cases\n",length(data3$b_submitdate)),
             paste0("At least 2 doses of Epinephrine\nn = ", length(which(data3$d_adre_no>=2))),
             paste0("Severe reaction \nrequiring hospitalisation\nn = ",length(which(refractory == "yes"))),
             paste0("Treatment with at least\ntwo 2nd line drugs\nn = ", length(which(refractory.strict=="yes"&refractory.ICU=="no"))),
             paste0("Final number of cases included in the study\nn = ",length(which(refractory.f=="yes"))-R1-R2-R3),
             paste0("Fatal reaction\nn = ",length(which(refractory.death=="yes"))),
             paste0("Reactions requiring\nICU\nn = ",length(which(refractory.ICU=="yes"&refractory.death=="no"))),
             paste0("Cases eliminated after manual revision:\nMultiple elicitors (n = 1)\n",
                    "Biphasic reactions responsive to Adrenaline (n = ", R1-1,")\n",
                    "Long time between adrenalin doses indicating responiveness (n = ",R2,")\n",
                    "Adequate response after second dose of adrenaline (n = ",R3,")")),
    #color = c("red", "green",
    #          "grey", "blue"),
    #value = c(3.5, 2.6, 9.4, 2.7),
    shape= "rectangle",
    width = c(1.8,3.5,1.8,1.8,3.5,1.8,1.8,5.8),
    x = c(1,1,1,3,1,-1, 1, 1),
    y = c(0,-1,-2,-3,-5,-3,-3,-4),
    height = c(rep(0.6,7),1))

edges <- create_edge_df(from = c(1,2,3,4,2,6,3,7,8),
                        to = c(2,3,4,8,6,8,7,8,5))
# Add the node data frame to the
# graph object to create a graph
# with nodes
create_graph() %>%
  add_node_df(node_df = nodes) %>%
  add_edge_df(edge_df = edges) %>%
  add_global_graph_attrs(
    attr = "splines",
    value = "ortho",
    attr_type = "graph") %>%
  add_global_graph_attrs(value="black",
                         attr = "fontcolor",
                         attr_type = "node") %>% #render_graph()
  export_graph(file_name = "analysis/figures/flow.png", file_type = "png", title = NULL,
               width = NULL, height = NULL)


