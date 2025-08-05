# Loading libraries
library(tidyverse)
library(mirt)
library(psych)
library(lessR)

## Loading the data
library(haven)
resilience <- read_sav("resilience.sav")


## selecting only USA##
USA.data <- resilience %>% 
  filter(CNTRYID==840) 


## deleting missing observations ##
USA.clean <- na.omit(USA.data)


## Special Function to Identify dif Items
get.dif.items <- function(f.data,p.val=.05,parms){
  r.warnings = ""
  keep.vars <- c("X2", "df", "p") 
  f.data <- f.data[keep.vars]
  f.data$p = round(f.data$p,3)
  if(missing(f.data)) return('Missing model output out.list')
  f.data$sig <- ifelse(f.data$p < p.val,'dif','no_dif')
  if(!missing(parms)){
    if(nrow(f.data) == nrow(parms)){
      f.data <- cbind(f.data,parms)
    }else{
      r.warnings = "There number of item parameters doesn't match the number of items "
      r.warnings = paste(r.warnings,"given to get.dif.items. Item parameters omitted.")
    }
  }
  dif.items <- subset(f.data, sig == 'dif')
  no.dif.items <- subset(f.data, sig == 'no_dif')
  if(!missing(parms) && nrow(f.data) == nrow(parms)){
    if(nrow(no.dif.items)>1){
      no.dif.items <- no.dif.items[order(-no.dif.items$a1),]
    }
  }
  r.list <- list(dif_items = dif.items, no_dif = no.dif.items, warnings = r.warnings)
  return(r.list)
}

## Identify Groups
## Dividing the data based on gender ##
Female <- USA.clean %>% filter(gender==1) %>% select (3:7)
Male <- USA.clean %>% filter(gender==2) %>% select (3:7)
FemaleN = 2268
MaleN <- 2280

## NB: I tried male as both "ref" and "foc" and got the same answer.
group <- c(rep('Ref', FemaleN), rep('Foc', MaleN))
dat <- rbind(Female, Male)

USA.clean %>% 
  group_by(gender, q4) %>%
  summarise(freq_q1 = n())

############# Step 1: Testing assumptions ############
####### Various Procedures to Check Dimensionality ######

## Parallel Analysis
fa.parallel(Female)
fa.parallel(Male)

## Unidimensionality Assessment
unidim(Female, cor = "poly")
unidim(Male, cor = "poly")

## Assessing GRM Models
## Reference Group
ref.model <- mirt(Female, model = 1, itemtype = "graded", SE=TRUE, verbose = FALSE)
M2(ref.model, type = "C2")
summary(ref.model)

ref.model.2 <- mirt(Female, model = 2, itemtype = "graded", SE=TRUE, verbose = FALSE)
M2(ref.model.2, type = "C2")
summary(ref.model.2)


## Anova test for the three models of the reference group
anova(ref.model, ref.model.2)

## Focal Group
foc.model <- mirt(Male, model = 1, itemtype = "graded", SE=TRUE, verbose = FALSE)
M2(foc.model, type = "C2")
summary(foc.model)

foc.model.2 <- mirt(Male, model = 2, itemtype = "graded", SE=TRUE, verbose = FALSE)
M2(foc.model.2, type = "C2")
summary(foc.model.2)

## Anova test for the three models of the focal group
anova(foc.model, foc.model.2)

## Assessing Item Fit
(Female.fit <- itemfit(ref.model))
(Male.fit <- itemfit(foc.model))

## Checking Frequencies
apply(Female, 2, table)
apply(Male, 2, table)

## Polychoric correlations
 polychoric(Female)
 polychoric(Male)
############# Step 2: Get Item Parameters ############
######### Freely estimated model ####################

model.free <- multipleGroup(dat, 1, group, verbose = FALSE)
coef(model.free, IRTpars = TRUE, printSE = TRUE, simplify = TRUE)  

############# Step 3: Baseline Model ##################

model.constrained <- multipleGroup(dat, 1, group, invariance = c(colnames(dat), 
                                                                 'free_means', 'free_var'))
(constrained.parameters <- coef(model.constrained, simplify = TRUE)[[1]][[1]]) 

############# Step 4: First round of DIF analyses (LRTs) - All Others As Anchors #################

(dif.drop <- DIF(model.constrained, c('a1','d1','d2','d3'), scheme = 'drop', 
                 seq_stat = .05, p.adjust = 'BH'))
get.dif.items(f.data=dif.drop,p.val=.05,parms=constrained.parameters)

############# Step 5: Specify a New Baseline Model Using Anchor Items (MaxA5 Approach)  #################

itemnames <- colnames(dat)
anc.items.names <- itemnames[c(3)]
test.items <- c(1,2,4,5)

model_anchor <- multipleGroup(dat, model = 1, group = group,
                              invariance = c(anc.items.names, 'free_means', 'free_var'))

(anchor.parms <-coef(model_anchor, simplify = TRUE)[[1]][[1]])


############# Step 6: Run first Final Invariance Tests  #################

(dif.anchor <- DIF(model_anchor, c('a1','d1','d2','d3'), items2test = test.items, p.adjust = 'BH'))
get.dif.items(f.data=dif.anchor,p.val=.05)



############# Step 7: Specify a New Baseline Model Using Anchor Items   #################
#### The previous results show that item 4 do not show dIF, #########
## Hence, specify a new baseline model with 3 and 4 as anchors ####

itemnames <- colnames(dat)
anc.items.names <- itemnames[c(3,4)]
test.items <- c(1,2,5)

model_anchor <- multipleGroup(dat, model = 1, group = group,
                              invariance = c(anc.items.names, 'free_means', 'free_var'))

(anchor.parms <-coef(model_anchor, simplify = TRUE)[[1]][[1]])


############# Step 6: Run second Final Invariance Tests  #################

(dif.anchor <- DIF(model_anchor, c('a1','d1','d2','d3'), items2test = test.items, p.adjust = 'BH'))
get.dif.items(f.data=dif.anchor,p.val=.05)


############# Step 9: Compute Effect Sizes  #################

## DTF
empirical_ES(model_anchor) 
empirical_ES(model_anchor, DIF=FALSE) 

##DIF
empirical_ES(model_anchor, DIF=FALSE, plot=TRUE) 
empirical_ES(model_anchor, plot=TRUE, as.table = TRUE) 

## DTF Plot (Gray Scale)
ES.item <- empirical_ES(model_anchor, plot=TRUE, as.table = TRUE)  
ES.item$main <- "" 
ES.item$par.settings$strip.background$col[[1]] <-"#E0E0E0" 
ES.item$panel.args.common$col[1:2] <- c('black', 'gray')  
ES.item 


## DIF Plots (Gray Scale)
ES.test <- empirical_ES(model_anchor, DIF=FALSE, plot=TRUE) 
ES.test$main <- "" 
ES.test$legend$top$args$key$text[[1]] <- c('Focal - Males', 'Reference - Females') 
ES.test$legend$top$args$key$points$col[1:2] <- c('black', 'gray') 
ES.test$panel.args.common$col[1:2] <- c('black', 'gray') 
ES.test 

## Article plots for items - tiff format
#tiff("Figure S3 - Expected Item Scores.tiff", width = 6, height = 4.5, units = 'in', res = 300)
#ES.item
#dev.off()

## Article plot for the scale
#tiff("Figure S4 - Expected Scale Scores.tiff", width = 6, height = 4.5, units = 'in', res = 300)
#ES.test
#dev.off()

