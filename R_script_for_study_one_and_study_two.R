## clean the environment -------------------------------------------------------
rm(list = ls(all=T)) # remove all variables (all=T takes care of hidden objects)
graphics.off() # turn off plots
cat("\014") # clear the console
## -----------------------------------------------------------------------------

################################################################################
#########################           NOTE          ##############################
################################################################################
################################################################################
###                                                                  ###########
###       This script only reproduces the core psychometric analysis ###########
###                                                                  ###########
###       The data has been striped and only essential analysis are  ###########
###          reported in this script. Therefore, information about   ###########
###        descriptive statistics and the like are not included to   ###########
###                maintain confidentiality with respondents         ###########
###                                                                  ###########
################################################################################
################################################################################


# Load necessary libraries -----------------------------------------------------
library(readxl)
library(psych)
library(ggfortify)
library(parameters)
library(grid)
library(directlabels)
library(gmodels)
library(DescTools)
library(misty)
library(lavaan) # for CFA
library(lavaanPlot) # for CFA
library(apaTables)
library(semTools)
library(tidyverse)

options(scipen = 999) #show as integers instead of exponential notations
options(digits=4) #Number of digits in output
# ------------------------------------------------------------------------------



## Data import and sanity checks -----------------------------------------------

d <- read_excel("study_one.xlsx")
summary(d)
str(d)

## -----------------------------------------------------------------------------


## Internal reliability estimates for auxiliary measures -----------------------

### GAD ###
# GAD screening is column 36:42
gad <- d[36:42]
omega(gad, plot = F)

# PHQ screening is column 43:51
phq <- d[43:51]
omega(phq, plot = F)


# RATHUS columns 6:35
# The items on the RAS that need to be reversed are: 
# 1, 2, 4, 5, 9, 11, 12, 13, 14, 15, 16, 17, 19, 23, 24, 26, and 30.
# However, they have already been reversed by means of the set-up of the online trial
rathus <- d[c(6:35)]

# rename rathus scale (RAS)
new_namesRAS <- paste0("ras", 1:30)
names(rathus)[1:30] <- new_namesRAS

omega(rathus, plot = F)


## QWB screening ---------------------------------------------------------------
# QWB is in colums 52:69

## swb = QWB in the code!
swb <- d[, c(52:69)]
summary(swb)

# quick look at omega and alpha

omega(swb[, c(1:18)], n.iter = 10, plot = F)
## -----------------------------------------------------------------------------


################# FACTOR ANALYSIS STUDY ONE ###### -----------------------------


swb %>% KMO() # .95 => suitable for FA

swb %>% bartlett.test() # significant, p<0.001 => data suitable for FA




# Make a scree plot ------------------------------------------------------------
scree <- fa.parallel(swb, 
                     fa="fa", 
                     main= "", 
                     n.iter=20,
                     quant=.95,
                     cor="cor",
                     fm = "ml",
                     plot=TRUE,
                     correct=.5, 
                     ylabel = "Eigenvalues", 
                     show.legend = F) 


# make the scree plot pretty ---------------------------------------------------
obs <- data.frame(scree$fa.values)
obs$type = c('Observed Data')
obs$num = c(row.names(obs))
obs$num = as.numeric(obs$num)
colnames(obs) = c('eigenvalue', 'type', 'num')

# Calculate quantiles for eigenvalues, but only store those from simulated CF model in percentile1
percentile = apply(scree$values,2,function(x) quantile(x,.95))
min = as.numeric(nrow(obs))
min = (4*min) - (min-1)
max = as.numeric(nrow(obs))
max = 4*max
percentile1 = percentile[min:max]

# Create data frame called sim with simulated eigenvalue data
sim = data.frame(percentile1)
sim$type = c('Simulated Data (95th %ile)')
sim$num = c(row.names(obs))
sim$num = as.numeric(sim$num)
colnames(sim) = c('eigenvalue', 'type', 'num')

#Merge the two data frames (obs and sim) together into data frame called eigendat
eigendat = rbind(obs,sim)


# APA theme for scree plot
apatheme=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text=element_text(family='Arial'),
        legend.title=element_blank(),
        legend.position=c(.7,.8),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'))

# Use data from eigendat. Map number of factors to x-axis, eigenvalue to y-axis, 
# and give different data point shapes depending on whether eigenvalue is observed or simulated
p = ggplot(eigendat, aes(x=num, y=eigenvalue, shape=type)) +
  #Add lines connecting data points
  geom_line()+
  #Add the data points.
  geom_point(size=4)+
  #Label the y-axis 'Eigenvalue'
  scale_y_continuous(name='Eigenvalue')+
  #Label the x-axis 'Factor Number', and ensure that it ranges from 1-max # of factors, increasing by one with each 'tick' mark.
  scale_x_continuous(name='Factor Number', breaks=min(eigendat$num):max(eigendat$num))+
  #Manually specify the different shapes to use for actual and simulated data, in this case, white and black circles.
  scale_shape_manual(values=c(16,1)) +
  #Add vertical line indicating parallel analysis suggested max # of factors to retain
  # geom_vline(xintercept = scree$nfact, linetype = 'dashed')+
  # kaiser criterion
  geom_hline(yintercept = 1, linetype = 'dashed', show.legend = T)+
  #Apply our apa-formatting theme
  apatheme

#Call the plot. 
p


## correlation between items ---------------------------------------------------

library(corrplot)
library(gtsummary)
# Correlation matrix of items
# sum(is.na(swb))
swb %>% tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd}) [{min}-{max}]"), 
                    digits = list(all_categorical() ~ 1))
660-166 # 494 == complete cases!


cormat <- swb %>% cor(use = "complete.obs")


## TABLE 4 ----------------------------------------
## Rationale for table 4:
# the difficulty parameter, one easy way to approximate it  is to compute the 
# mean of each item, and rankorder them and correlate them to see if we have a 
# level discrimination problem with this solution.


# The swb data is stored in a data frame called swb
swb_means <- swb %>%
  summarise(across(starts_with("swb"), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = starts_with("swb"), names_to = "item", values_to = "mean")

# Rank ordering based on means
swb_means <- swb_means %>%
  mutate(rank = rank(-mean)) # Negative for descending order (easiest to hardest)

swb_means
# Correlating rank orders with item means
correlation <- cor(swb_means$rank, swb_means$mean)

cor.test(swb_means$rank, swb_means$mean)
# Output the correlation
print(correlation)


### ----------------------------------------------------------------------------


## EFA study one ---------------------------------------------------------------
efa1 <- fa(swb, 
           nfactors = 1, 
           rotate="oblimin", 
           fm="ml", 
           missing=TRUE, 
           impute="mean", ) 
print(efa1, sort = T, cut = 0.3)


# function for creating FA table

fa_table <- function(x, cut) {
  #get sorted loadings
  loadings <- fa.sort(x)$loadings %>% round(3)
  #supress loadings
  loadings[loadings < cut] <- ""
  #get additional info
  add_info <- cbind(x$communalities, 
                    x$uniquenesses,
                    x$complexity) %>%
    # make it a data frame
    as.data.frame() %>%
    # column names
    rename("Communality" = V1,
           "Uniqueness" = V2,
           "Complexity" = V3) %>%
    #get the item names from the vector
    rownames_to_column("item")
  #build table
  loadings %>%
    unclass() %>%
    as.data.frame() %>%
    rownames_to_column("item") %>%
    left_join(add_info) %>%
    mutate(across(where(is.numeric), round, 2))
}

fa_table(efa1, 0.3)

efa1$values


alpha(swb)
omega(swb, plot = F, nfactors = 1, n.iter=10)


# convergent validity ----------------------------------------------------------

library(lavaan)
library(lavaanPlot)
library(semTools)
library(semPlot)
library(psych)
library(readxl)

X <- cbind.data.frame(swb, rathus)

Fakt_swb_ras <- 'swb =~ swb1+swb2+swb3+swb4+swb5+swb6+swb7+swb8+swb9+swb10+swb11+
                        swb12+swb13+swb14+swb15+swb16+swb17+swb18
                 ras =~ ras1+ras2+ras3+ras4+ras5+ras6+ras7+ras8+ras9+ras10+ras11+
                        ras12+ras13+ras14+ras15+ras16+ras17+ras18+ras19+ras20+ras21+
                        ras22+ras23+ras24+ras25+ras26+ras27+ras28+ras29+ras30'

fit <- sem(Fakt_swb_ras, data = X)

summary(fit, standardized=TRUE)

fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea","TLI","aic","bic"))

semPaths(fit,  "std", )

# Figure S1: 
semPlot::semPaths(fit, "std", layout = "tree2", sizeMan = 2)
semPlot::semPaths(fit, 
                  "std", 
                  layout = "tree2",
                  # nCharNodes =  10,
                  intercepts = T,
                  residuals = T,
                  nDigits = 1,
                  edge.color = "black",
                  color = "lightgray",
                  what = "path",
                  rotation = 1,
                  sizeMan = 2,
                  label.cex = 1.5,  # Adjust this value as needed for larger font size
                  nodeLabels = c(1:18, 1:30,
                                 # "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                 # "11", "12", "13", "14", "15", "16", "17", "18",
                                 # "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                 # "11", "12", "13", "14", "15", "16", "17", "18",
                                 # "19", "20", "21", "22", "23", "24", "25", "26", "27", "28",
                                 #  "27", "28", "29", "30",
                                 "QWB", "RAS"),
                  mar = c(5, 1, 5, 1)
)

reliability(fit)

discriminantValidity(fit)

discriminantValidity(fit, merge = TRUE)



## Criterion-Related Validity (study one) --------------------------------------

# make sum scores numeric variables
d$RATHUS_sceening_SUM <- d$RATHUS_sceening_SUM %>% as.numeric()
d$QWB_screening_SUM <- d$QWB_screening_SUM %>% as.numeric()
d$`GAD-7_screening_SUM` <- d$`GAD-7_screening_SUM` %>% as.numeric()
d$`PHQ-9_screening_SUM` <- d$`PHQ-9_screening_SUM` %>% as.numeric()

# QWB and RAS scale:
cor.test(d$RATHUS_sceening_SUM, d$QWB_screening_SUM)

# QWB and GAD-7
cor.test(d$`GAD-7_screening_SUM`, d$QWB_screening_SUM)

# QWB and PHQ-9
cor.test(d$`PHQ-9_screening_SUM`, d$QWB_screening_SUM)

### ----------------------------------------------------------------------------


## retest analysis ------------------------------------------------------------

### Combine datasets to calculate ICC3 for test-retest reliability 

# halfway points
d$QWB_SUM_halftime <- d$QWB_SUM_halftime %>% as.numeric()
sbwRetestM <- cbind.data.frame(d$QWB_screening_SUM, d$QWB_SUM_halftime)

# post-assessment
d$QWB_SUM_post <- d$QWB_SUM_post %>% as.numeric()
sbwRetestP <- cbind.data.frame(d$QWB_screening_SUM, d$QWB_SUM_post)

# one year follow-up
d$QWB_SUM_one_year_followup <- d$QWB_SUM_one_year_followup %>% as.numeric()
sbwRetestY <- cbind.data.frame(d$QWB_screening_SUM, d$QWB_SUM_one_year_followup)



## ICC3 proper 

# ICC for pre test and middle of intervention measurement
psych::ICC(x = sbwRetestM, missing = T, alpha = .05, lmer = T) # ICC3 = 0.65 [95% CI:  0.61, 0.69]

# ICC for pre and post measurement
psych::ICC(x = sbwRetestP, missing = T, alpha = .05, lmer = T) # ICC3 = 0.50 [95% CI:  0.45, 0.54]

# ICC for pre test and one-year follow up 
psych::ICC(x = sbwRetestY, missing = T, alpha = .05, lmer = T) # ICC3 = 0.51 [95% CI:  0.46, 0.55]


# look at post test and one-year follow up (not reported in paper)
sbwRetestPY <- cbind.data.frame(d$QWB_SUM_post, d$QWB_SUM_one_year_followup)
psych::ICC(x = sbwRetestPY, missing = T, alpha = .05, lmer = T) # ICC3 = 0.51 [95% CI:  0.46, 0.55]

### ----------------------------------------------------------------------------


## Cutoff point analysis -------------------------------------------------------


# Load necessary packages
library(pROC)

# Create a binary outcome based on PHQ-9 and GAD-7
binary_outcome <- ifelse(d$`PHQ-9_screening_SUM` >= 10 | d$`GAD-7_screening_SUM` >= 8, 1, 0)

# Perform ROC analysis
roc_obj <- roc(binary_outcome, d$QWB_screening_SUM, 
               ci=T,
               # arguments for plot
               plot=TRUE, auc.polygon=F, max.auc.polygon=F, 
               # grid=TRUE,
               print.auc=TRUE)

# Get the coordinates for the 'best' cutoff point
cutoff_coords <- coords(roc_obj, "best", ret=c("threshold", "sensitivity", "specificity"), transpose = FALSE)

# Extract the cutoff value from the coordinates
cutoff_value <- cutoff_coords["threshold"]

# Add the point to the plot
points(cutoff_coords["specificity"], cutoff_coords["sensitivity"], col="red", pch=19, cex=1.5)

# Add a label for the point
text(cutoff_coords["specificity"], cutoff_coords["sensitivity"], 
     labels=sprintf("QWB cutoff: %s", cutoff_value), pos=2.5)
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------


################################################################################
####################                 ###########################################
####################    STUDY TWO    ###########################################
####################                 ###########################################
################################################################################




## Study two ------------------------------------------------------------------


## clean the environment -------------------------------------------------------
rm(list = ls(all=T)) # remove all variables (all=T takes care of hidden objects)
graphics.off() # turn off plots
cat("\014") # clear the console
## -----------------------------------------------------------------------------


# Read in study two data -------------------------------------------------------
dd <- read_excel("study_two.xlsx")
summary(dd)
str(dd)

## -----------------------------------------------------------------------------

               #########################################
               #####                         ###########
               #####     NOTE: swb = QWB     ###########
               #####                         ###########
               #########################################

# swb items extract
swbNY <- dd[, c(2:19)]



### ----------------------------------------------------------------------------


#### SWB sum score ------------------

dd$swb_SUM <- dd$swb1 + dd$swb2 + dd$swb3 + dd$swb4 + dd$swb5 + dd$swb6 +
  dd$swb7 + dd$swb8 + dd$swb9 + dd$swb10 + dd$swb11 + dd$swb12 +
  dd$swb13 + dd$swb14 + dd$swb15 + dd$swb16 + dd$swb17 + dd$swb18
sum(is.na(dd$swb_SUM))

# Internal reliability for QWB -------------------------------------------------
omega(swbNY, nfactors = 1, plot = F, n.iter = 10)


### ----------------------------------------------------------------------------

## bbq sum score ---------------------------------------------------------------

# Item 1 * Item 2 = leisure_time
dd$leisure_time <- dd$bbq1 * dd$bbq2

# Item 3 * Item 4 = view_on_life
dd$view_on_life <- dd$bbq3 * dd$bbq4

# Item 5 * Item 6 = creativity
dd$creativity <- dd$bbq5 * dd$bbq6

# Item 7 * Item 8 = learning
dd$learning <- dd$bbq7 * dd$bbq8

# Item 9 * Item 10 = friends_and_friendship
dd$friends_and_friendship <- dd$bbq9 * dd$bbq10

# Item 11 * Item 12 = view_of_self
dd$view_of_self <- dd$bbq11 * dd$bbq12

dd$BBQ_SUM <- dd$leisure_time + dd$view_on_life +
  dd$creativity + dd$learning +
  dd$friends_and_friendship + dd$view_of_self

### ----------------------------------------------------------------------------


## Sanity check ----------------------------------------------------------------
# sanity check for scree plot in study two (not reported in the paper)
fa.parallel(swbNY, 
            fa="fa", 
            main= "", 
            n.iter=20,
            quant=.95,
            cor="cor",
            fm = "ml",
            plot=TRUE,
            correct=.5, 
            ylabel = "Eigenvalues", 
            show.legend = F) # suggests 6 factors


### CFA ------------------------------------------------------------------------

## fit all data to one factor

onefactor <- 'f1 =~ swb1 + swb2 + swb3 + swb4 + swb5 + swb6 + swb7 + swb8 +
                    swb9 + swb10 + swb11 + swb12 + swb13 + swb14 + swb15 + 
                    swb16 + swb17 + swb18'

# Fit the model to the data
cfamodel <- sem(model = onefactor, data = dd, estimator = "WLSMV")
cfamodel%>% summary(standardized=T, ci=F, fit.measures=TRUE, )

standardizedsolution(cfamodel)
inspect(cfamodel, 'r2')

# figure 3:
semPlot::semPaths(cfamodel, 
                  "std", 
                  layout = "tree2",
                  nCharNodes =  10,
                  intercepts = T,
                  residuals = T,
                  nDigits = 2,
                  edge.color = "black",
                  color = "lightgray",
                  what = "path",
                  rotation = 2,
)

# modification indicies 
modindices(cfamodel, sort. = T)

## evaluation of the mod indecies (not reported in the paper)

# onefactor_modified <- '
#   f1 =~ swb1 + swb2 + swb3 + swb4 + swb5 + swb6 + swb7 + swb8 +
#         swb9 + swb10 + swb11 + swb12 + swb13 + swb14 + swb15 +
#         swb16 + swb17 + swb18
#   swb4 ~~ swb5
#   swb5 ~~ swb12
#   swb1 ~~ swb2
#   swb2 ~~ swb16
#   swb12 ~~ swb13
#   swb11 ~~ swb17



# Polychoric matrix ------------------------------------------------------------
poly_gogn <- polychoric(dd[, c(3:20)], smooth = T)
rho <- poly_gogn$rho

# for colors 
gr <- colorRampPalette(c("#B52127", "white", "#2171B5"))


cor.plot(poly_gogn$rho, 
         numbers=T, 
         upper=FALSE, 
         main = "Polychoric Correlation Matrix", 
         show.legend = FALSE, 
         cex = 0.5, 
         # colors = gr
         gr = gr,labels = c(1:18))





## -----------------------------------------------------------------------------¨¨


## cireterion related validity in study two ------------------------------------


# criterion-related validity BBQ and SWB ---------------------------------------
cor.test(dd$swb_SUM, dd$BBQ_SUM)



library(lavaan)
library(lavaanPlot)
library(semTools)
library(semPlot)
library(psych)
library(readxl)

swb <- dd[, 2:19]
bbq <- dd[, 33:38]

X <- cbind.data.frame(swb, bbq)

Fakt_swb_bbq <- 'swb =~ swb1+swb2+swb3+swb4+swb5+swb6+swb7+swb8+swb9+swb10+swb11+
                        swb12+swb13+swb14+swb15+swb16+swb17+swb18
                 bbq =~ leisure_time + view_on_life + creativity + learning + 
                        friends_and_friendship + view_of_self'

fit <- sem(Fakt_swb_bbq, data = X)

# sanity checks 
summary(fit, standardized=TRUE)
fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea","TLI","aic","bic"))


# figure s2:
semPlot::semPaths(fit, 
                  "std", 
                  layout = "tree2",
                  # nCharNodes =  10,
                  intercepts = T,
                  residuals = T,
                  nDigits = 1,
                  edge.color = "black",
                  # edge.label.cex = 0.8, 
                  sizeMan = 3,  
                  color = "lightgray",
                  what = "path",
                  rotation = 1, 
                  mar = c(5, 5, 5, 3), 
                  nodeLabels = c("1", "2", "3", "4", "5", "6",
                                 "7", "8", "9", "10","11", "12",
                                 "13", "14", "15", "16","17", "18",
                                 "LT", "VoL", "C", "L", "FaF", "VoS", "QWB", "BBQ"),
                  label.cex = 1.5  
)

reliability(fit)

discriminantValidity(fit)

discriminantValidity(fit, merge = TRUE)





# internal reliability of BBQ in study two --------
omega(bbq, plot = F)


