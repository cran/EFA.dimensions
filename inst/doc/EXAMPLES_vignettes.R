## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "",
  fig.width = 7,
  fig.align = "center",
  background='#FFFFFF'
)

## -----------------------------------------------------------------------------
library(EFA.dimensions)

## ----collapse=FALSE-----------------------------------------------------------
ESEM(data = data_HS_1939[,-10], method = 'efa_blocks', Nfactors = 3)

## ----collapse=FALSE-----------------------------------------------------------
ESEM(data = data_HS_1939[,-10], method = 'startvalues', Nfactors = 3)

## ----collapse=FALSE-----------------------------------------------------------
ESEM(data = data_HS_1939[,-10], method = 'efa_blocks', target_keys = c(1,1,1,2,2,2,3,3,3))

## ----collapse=FALSE-----------------------------------------------------------
ESEM(data = data_HS_1939[,-10], method = 'startvalues', target_keys = c(1,1,1,2,2,2,3,3,3))

## ----collapse=FALSE-----------------------------------------------------------
ESEM(data = data_SDT, method = 'efa_blocks', target_keys = c(1,1,1,1,2,2,2,2,3,3,3,3))

## ----collapse=FALSE-----------------------------------------------------------
MPlus_4.1b_data <- read.table("https://www.statmodel.com/usersguide/chap4/ex4.1b.dat", header=FALSE)
colnames(MPlus_4.1b_data) <- paste0("y", 1:ncol(MPlus_4.1b_data))

ESEM(data = MPlus_4.1b_data, method = 'efa_blocks', Nfactors = 4, estimator = 'ML')

# compare with the MPlus output at  
# https://www.statmodel.com/usersguide/chap4/ex4.1part2.html

## ----collapse=FALSE-----------------------------------------------------------
BIFACTOR(rawdata = data_HS_1939[,-10],
         bifactor_kind = 'ESEM',
         Nfactors = 4,
         LV_options = list(group_keys = NULL,
                           estimator = 'MLR',
                           rotation = 'bigeomin',
                           resid_correls = NULL,
                           LV_names = NULL,
                           ordered = FALSE))


## ----collapse=FALSE-----------------------------------------------------------
BIFACTOR(rawdata = data_SDT,   
         bifactor_kind = 'CFA', 
         LV_options = list(group_keys = c(1,1,1,1,2,2,2,2,3,3,3,3), 
                           estimator='MLR', 
                           resid_correls=NULL, 
                           LV_names=NULL,
                           ordered = FALSE))

## ----collapse=FALSE-----------------------------------------------------------
BIFACTOR(rawdata = data_SDT,   
         bifactor_kind = 'ESEM', 
         EFA_options   = list(extraction = 'minres', rotation = 'oblimin', Nfactors = 3), 
         LV_options = list(group_keys = c(1,1,1,1,2,2,2,2,3,3,3,3), 
                           estimator='MLR', 
                           resid_correls=NULL, 
                           LV_names=NULL,
                           ordered = FALSE))

## ----collapse=FALSE-----------------------------------------------------------
MPlus_4.7_data <- read.table("https://www.statmodel.com/usersguide/chap4/ex4.7.dat", header=FALSE)
colnames(MPlus_4.7_data) <- paste0("y", 1:ncol(MPlus_4.7_data))

BIFACTOR(rawdata = MPlus_4.7_data,   
         Nfactors = 3,
         bifactor_kind = 'ESEM', 
         LV_options = list(group_keys = NULL,
                           estimator='ML',
                           rotation = 'bigeomin',
                           resid_correls=NULL,
                           LV_names=NULL,
                           ordered = FALSE))

# compare with the MPlus output at  
# hhttps://www.statmodel.com/usersguide/chap4/ex4.7.html

## ----collapse=FALSE-----------------------------------------------------------
HS_1939_model <- 'visual  =~ test1 + test2 + test3
                  textual =~ test4 + test5 + test6
                  speed   =~ test7 + test8 + test9'

HS_1939_output <-
  Factorial_Invariance(data = data_HS_1939, group = 'school', LV_model = HS_1939_model)

PLOT_Invariance(model_object = HS_1939_output,
                invar_model = 'Scalar',
                plot_types = c('loadings', 'LV_distribs')) 


# for more detailed statistical comparisons of the latent variable means, try
# the GROUP.DIFFS function from the DFA.CANCOR package, as follows:
# install.packages(DFA.CANCOR); library(DFA.CANCOR)
# 
# LV_scores <- HS_1939_output$LV_scores
# 
# GROUP.DIFFS(data = LV_scores, GROUPS = 'school', DV = 'visual')
# 
# GROUP.DIFFS(data = LV_scores, GROUPS = 'school', DV = 'speed')

## ----collapse=FALSE-----------------------------------------------------------
HS_1939_output <-
Factorial_Invariance(data = data_HS_1939, group = 'school',
                     LV_keys = c(1,1,1, 2,2,2, 3,3,3),
                     LV_names = c('visual', 'textual', 'speed'))

PLOT_Invariance(model_object = HS_1939_output, 
                invar_model = 'Configural',
                plot_types = c('loadings', 'ints_slopes')) 

## ----collapse=FALSE-----------------------------------------------------------
Factorial_Invariance(data = data_HS_1939, group = 'school',
                     LV_keys = c(test1 = 1, test2 = 1, test3 = 1,
                                 test4 = 2, test5 = 2, test6 = 2,
                                 test7 = 3, test8 = 3, test9 = 3) )

## ----collapse=FALSE-----------------------------------------------------------
# 2015 Brown - Confirmatory Factor Analysis for Applied Research  p 247, Table 7.11

# Tests of measurement invariance and population heterogeneity of DSM-IV
# major depressive disorder in men and women
data_Brown_2015 <- read.table("http://people.bu.edu/tabrown/Ch7/MDDALL.dat")
names(data_Brown_2015) <- c("sex", paste("mdd", 1:9, sep = ""))
data_Brown_2015$sex <- factor(data_Brown_2015$sex, levels = c(0, 1), labels = c("female", "male"))

# using LV_model
model.mdd <- '
MDD =~ mdd1 + mdd2 + mdd3 + mdd4 + mdd5 + mdd6 + mdd7 + mdd8 + mdd9
mdd1 ~~ mdd2
'
Factorial_Invariance(data = data_Brown_2015, group = 'sex', LV_model = model.mdd)

# using LV_keys
Brown_2015_output <-
  Factorial_Invariance(data = data_Brown_2015, group = 'sex', 
                     LV_keys = c(1,1,1,1,1,1,1,1,1),
                     LV_resid_correls = c('mdd1 ~~ mdd2'))

PLOT_Invariance(model_object = Brown_2015_output,
                invar_model = 'Scalar',
                plot_types = c('loadings', 'LV_distribs')) 

