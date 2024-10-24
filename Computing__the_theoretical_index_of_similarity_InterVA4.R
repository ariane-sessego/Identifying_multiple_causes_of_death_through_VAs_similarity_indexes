##### Computing the theoretical index of similarity based on InterVA-4 ####

#### Script for the replicability of the results presented in
### "Studying multiple causes of death in the absence of death certificates:
### taking advantage of verbal autopsies"

#### Ariane SESSEGO under the supervision of Geraldine DUTHE (INED, France), 
#### in collaboration with Bruno Lankoandé and Dianou Kassoum (ISSP, Ouagadougou) 
#### March 2023 #####
#### R.version : 4.1.2

### This script is divided in 4 main sections :
### 1. The preparation of InterVA's probability matrix (probbase) for the computation of the indexes
### 2. The computation of the indexes and formating into matrixes and edgelists to be used for data analysis
### 3. Some simple visualisations of the results, and comparison of the indexes
### 4. Robustness checks with age-sex specific indexes and malaria and HIV prevalence specific indexes



### Preparation ----- 

## Packages 
library(questionr)
library(tidyverse)
library(haven)
library(dplyr)
library(tidyr)
library(foreign)
library(igraph)
library(netdiffuseR) 


# Set the working directory to file of the script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

#### 1. The preparation of InterVA's probability matrix (probbase) -----

## InterVA-4 requires to specify the levels of malaria and HIV/AIDs in the population before computing, 
## To adjust the a priori prevalence of these diseases. 
## We use here the the probability matrix (probbase) on the settings malaria "Very Low" (VL) and HIV/AIDs "Very Low" (VL)
## as they are the default settings of the model, but it can be computed for any specifications.

## It only necessitates a minor change in the a priori matrix, according to the setting chosen :
##  The a priori prevalence of the disease must be changed to the apropriate letter according to the setting chosen,
## The a priori prevalence of sickle cell disease is set the same value of the one for malaria.
## The settings translate in letters as follows High = B, Low = C and Very Low = E.

## These minor changes do not significantly affect the results of the index on the general level, 
### (as show in Sessego, Ariane, Géraldine Duthé, Bruno Lankoandé, et Dianou Kassoum. 2021. 
## « Studing Multiple Causes of Death in the Absence of Death Certificates : 
## Taking Advantage of Probabilistic Cause of Death Estimation Methods (InterVA-4) ». 
## Ined, Documents de travail, 268), so stuck with the probbase VL/VL.
## They might however be more important for indexes based on individual symptoms,
## the probbase would hence have to be changed accordingly.

### a. Loading  the probbase ----

prob <- read.csv("Input_probbase/probbaseVLVL_letters.csv")
prob$X <- NULL

### b. Translating the letters in numbers ----
# Transforming all the columns in character variables 
for (i in 2:61) { prob[,i] <- as.character(prob[,i])}
str(prob)
#### Transforming the letters into numbers (as presented in Byass, 2014) 
prob[prob == "I"] <- "1"
prob[prob == "A+"] <- "0.8"
prob[prob == "A"] <- "0.5"
prob[prob == "A-"] <- "0.2"
prob[prob == "B+"] <- "0.1"
prob[prob == "B"] <- "0.05"
prob[prob == "B-"] <- "0.02"
prob[prob == "B -"] <- "0.02"
prob[prob == "C+"] <- "0.01"
prob[prob == "C"] <- "0.005"
prob[prob == "C-"] <- "0.002"
prob[prob == "D+"] <- "0.001"
prob[prob == "D"] <- "0.0005"
prob[prob == "D-"] <- "0.0001"
prob[prob == "E"] <- "0.00001"
prob[prob == "N"] <- "0.0"

### c. Exporting the resulting database ----
 # write.csv(prob, "Input_probbase/probbaseVLVL_numbers.csv")


### 2. Computing the indexes ----
## a. creating the necessary functions -----

# Let A and B be vectors that associate the probability of presenting each symptom given cause a and b respectively

## The Euclidean index
indice_euclidien <- function(A,B){
  A <- as.matrix(as.numeric(A)) # Transforming the vectors into numeric vectors to allox for computation  
  B <- as.matrix(as.numeric(B))
  indice <- norm(A - B, type = "F")/ norm(A + B, type = "F") # type = "F" i.e. chosing the Euclidean norm
  return(indice)
}

# The absolute norm index
indice_norme1 <- function(A,B){
  A <- as.matrix(as.numeric(A)) # Transforming the vectors into numeric vectors to allow for computation  
  B <- as.matrix(as.numeric(B))
  indice <- norm(A - B, type = "O")/ norm(A + B, type = "O") # type = "0" i.e. chosing norm 1
  return(indice)
}

### Testing
# 1 maximum dissimilarity
indice_euclidien(c(0.01,1,0), c(0.1,1,0))
indice_norme1(c(0.01,1,0), c(0.1,1,0))

indice_euclidien(c(0,1,0), c(0,1,0))
# indice_norme1(c(0,1,0), c(0,1,0))
# indice_norme1(c(1,0,1), c(0,1,0))
# indice_euclidien(c(1,0,1), c(0,1,0))
# indice_norme1(c(1,1,0), c(0,1,0))
# indice_euclidien(c(1,1,0), c(0,1,0))
# 
# indice_norme1(prob$B_SEPSIS.C.7, prob$B_SEPSIS.C.7)
# indice_norme1(prob$B_SEPSIS.C.7, prob$B_PNEUM.C.7)


# Normalising on the number of indicators -----
indice_unif_normalisation <- function(A,B){
  A <- as.matrix(as.numeric(A)) # Transforming the vectors into numeric vectors to allow for computation  
  B <- as.matrix(as.numeric(B))
  indice <- norm(A - B, type = "F")/ length(A) #Number of indicators
  return(indice)
}
#1 == parfaitement differents
# Testing 
# indice_unif_normalisation(c(0,1,0), c(1,0,1))
# indice_unif_normalisation(c(1,0,1), c(1,0,1)) # Parfaitement égaux - indice de dissimilarité
# indice_unif_normalisation(c(1,1,1), c(1,0,1))
# indice_unif_normalisation(c(1,1,1), c(1,0.5,1))
# 


# Scalar product ----
indice_scalar_product <- function(A,B){
  A <- as.matrix(as.numeric(A)) # Transforming the vectors into numeric vectors to allow for computation  
  B <- as.matrix(as.numeric(B))
  indice <- t(A)%*%B/(norm(A, type = "F")*norm(B, type = "F") )
  return(1- indice[1,1]) # 1- the index, to have the 1 == maximally different, 0 perfectly similar, as for the others.
}
#Tests
# indice_scalar_product(c(0,1,0), c(1,0,1))
# 1-indice_scalar_product(c(1,0,1), c(1,0,1)) #
# 1-indice_scalar_product(c(1,1,1), c(1,0.5,1))





### b. systematising the computing for all possible combinations  -----

# Given a matrix M of probabilities associating to each cause (in the colums),
# probabilities to present each symptom (in rows), 
# this function allows to compute the index for every possible combination of two causes

matricer <- function(M, f){ ## M the matrix and f the index function of choice
  
  dimM<- as.matrix(dim(M))[2,1] ### store the number of columns i.e. possible causes
  
  comb <- combn(c(1:dimM), 2)# we compute all the possible combinations of columns
  dimcomb <- as.matrix(dim(comb))[2,1] # we put it in the right format 
  
  
  # We create an empty index matrix of the right size
  I <- matrix(data = 0, nrow = dimM, ncol = dimM)
  # We fill the matrix
  for (i in (1:dimcomb)){
    I[comb[1,i],comb[2,i]] = f(unlist(M[,comb[1,i]]),unlist(M[,comb[2,i]]))
  }
  # rename the colums 
  colnames(I) <- colnames(M)
  rownames(I) <- colnames(M)
  
  return(I) # return our index matrix !!
}

### c. computing the indexes ----

# we suppress the first colum as it is only the label of each indicator 
prob$INDIC.C.10 <- NULL

# we compute the matrix indexes of all the indicators
 
#  for (i in c("indice_euclidien",
#              "indice_norme1",
#              "indice_unif_normalisation",
#             "indice_scalar_product")){
# return(assign(paste0("M_", i),as.data.frame(matricer(prob, noquote(i)))))
#  }
# 
#  as.character(i)
# 
# i <- 3

  
# we compute the Euclidean index matrix
M_indice_euclidien <- as.data.frame(matricer(prob, indice_euclidien))

# we compute the absolute norm index matrix
M_indice_norme1 <- as.data.frame(matricer(prob, indice_norme1))

# indice_unif_normalisation
M_indice_unif_normalisation <- as.data.frame(matricer(prob, indice_unif_normalisation))

# we compute the absolute norm index matrix
M_indice_scalar_product <- as.data.frame(matricer(prob, indice_scalar_product))


## Exporting these matrixes
# write.csv(M_indice_euclidien, "Output_indexes/VLVL_Matrix_Eucl_index.csv")
# write.csv(M_indice_norme1, "Output_indexes/VVLVL_Matrix_norm1_index.csv")
# write.csv(M_indice_unif_normalisation, "Output_indexes/VVLVL_Matrix_normalise_index.csv")
# write.csv(M_indice_scalar_product, "Output_indexes/VVLVL_Matrix_scalar_product_index.csv")

### d. creating an edgelist database of the index ----

# The matrix format is practical for computation, but also for simple visualisation of the indexes.
# However, to then merge the index with VA data, an edgelist format 
# (with in colums for the association of causes and the associated index)
# is much more practical, hence this part of the script

# for this we use packages designed for Network analysis
library(igraph)
library(netdiffuseR) 

# However, to allow for good formating, especially for cause names, we will want to replace cause numbers by names.
### To do so we use the following table presenting the equivalences of designation in english and french
### very practical to translate tables and graphics.
### It was constructed from the causetext table already available in the InterVA-4 package
causetext <-  read.csv("Input_probbase/causetext_AS.csv", sep = ",")

# We suppress the unnecessary rows 
causetext <- causetext[4:63,]

# We create a function that allows us to add a colum with the equivalent designation according to our .xlsx table
renommer <- function(
    base_var, # data base where we want to rename a variable
    base_noms # base of the equivalent designations
    , arenommer_entreguill # name of the variable to be renamed (in "") 
    , equivalent # name of this designation in causetext
    , nv_nom # the name of the variable of the new desired designation
    , nvnom_guill # name of the new variable created
){
  equivalent <- as.double(equivalent)
  
  # we extract the two designations that we need and bind them together
  C <- bind_cols(equivalent, nv_nom)
  
  # we rename the columns 
  colnames(C) <- bind_cols(arenommer_entreguill, nvnom_guill)
  
  # we join the two tables
  base_var <- left_join(base_var, C, by = arenommer_entreguill )
  
  return(base_var)
}

## For Euclidean index

edgelist_euclidien <- adjmat_to_edgelist( # the very useful function from netdiffuseR, 
  # allowing to transform an adjacency matrix to an edgelist
  matricer(prob, indice_euclidien),
  undirected =FALSE,
  keep.isolates = TRUE
) %>% as.data.frame() %>% 
  renommer(causetext, "ego", causetext$numero_cause
           , causetext$CAUSETXT.C.40, "ego_name") %>% 
  renommer(causetext, "alter", causetext$numero_cause
           , causetext$CAUSETXT.C.40, "alter_name")

# Norm 1

edgelist_norme1 <- adjmat_to_edgelist(
  matricer(prob, indice_norme1),
  undirected =FALSE,
  keep.isolates = TRUE
) %>% as.data.frame() %>% 
  renommer(causetext, "ego", causetext$numero_cause
           , causetext$CAUSETXT.C.40, "ego_name") %>% 
  renommer(causetext, "alter", causetext$numero_cause
           , causetext$CAUSETXT.C.40, "alter_name")


# Normalised 

edgelist_normalised <- adjmat_to_edgelist(
  matricer(prob, indice_unif_normalisation),
  undirected =FALSE,
  keep.isolates = TRUE
) %>% as.data.frame() %>% 
  renommer(causetext, "ego", causetext$numero_cause
           , causetext$CAUSETXT.C.40, "ego_name") %>% 
  renommer(causetext, "alter", causetext$numero_cause
           , causetext$CAUSETXT.C.40, "alter_name")



# Scalar product


edgelist_scalar_product <- adjmat_to_edgelist(
  matricer(prob, indice_scalar_product),
  undirected =FALSE,
  keep.isolates = TRUE
) %>% as.data.frame() %>% 
  renommer(causetext, "ego", causetext$numero_cause
           , causetext$CAUSETXT.C.40, "ego_name") %>% 
  renommer(causetext, "alter", causetext$numero_cause
           , causetext$CAUSETXT.C.40, "alter_name")


###### We join the two in one single database 
edgelist_indices <- left_join(edgelist_norme1, edgelist_euclidien,
                               by = c("ego", "alter", "time", "ego_name", "alter_name") ) %>%
  rename("indice_norm" = "value.x") %>% 
  rename("indice_euclidien" = "value.y") %>% 
  left_join(edgelist_scalar_product,
            by = c("ego", "alter", "time", "ego_name", "alter_name") ) %>% 
  rename("indice_scalar_product" = "value") %>% 
  left_join(edgelist_normalised,
            by = c("ego", "alter", "time", "ego_name", "alter_name") ) %>% 
  rename("indice_normalised" = "value")



q_eucl <-quantile(edgelist_indices$indice_euclidien, probs = 0.25)
  library(dplyr)
library(stats)
edgelist_indices<- edgelist_indices %>% 
  mutate(coocc_eucl = case_when(indice_euclidien >= quantile(edgelist_indices$indice_euclidien, probs = 0.25) ~ "yes",
                                 T ~ "no"), 
         oocc_norm = case_when(indice_norm >= quantile(edgelist_indices$indice_norm, probs = 0.25) ~ "yes",
                               T ~ "no"), 
         oocc_scalar = case_when(indice_scalar_product >= quantile(edgelist_indices$indice_scalar_product, probs = 0.25) ~ "yes",
                               T ~ "no"), 
         oocc_normalised = case_when(indice_normalised >= quantile(edgelist_indices$indice_normalised, probs = 0.25) ~ "yes",
                               T ~ "no")
    )


# write.csv(edgelist_indices, "Output_indexes/VLVL_Edgelist_indexes.csv")






### 3. Simple visualisations -----


## Tidying the data frame ----

edgelist_indices_tidy <- pivot_longer(edgelist_indices,
             cols = c("indice_norm", 
             "indice_euclidien", 
             "indice_scalar_product", 
             "indice_normalised"),
             names_to = "type")

## a. Distribution of the indexes -----
library(RColorBrewer)
ggplot(edgelist_indices_tidy) +
  geom_freqpoly(mapping = aes(x = value, color = type))+
  scale_x_continuous(limits = c(0,1)) +
  labs(x = "Value of the index",
       y = "Frequency",
       title ="Distribution of the similarity indexes (InterVA4) ", 
       subtitle = "For all possible 1770 combination of causes",
       caption = "Computed from the default probability matrix of InterVA4 (probbase Malaria : VL, HIV : VL)"
  ) +
   scale_color_manual(name = "Type of index",
                      labels = c("indice_norm" = "Absolute norm",
                                 "indice_euclidien" = "Euclidean",
                                 "indice_scalar_product" = "Scalar product",
                                 "indice_normalised" = "Uniform normalisation"),
                      values = brewer.pal(4, "Set1")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))  
ggsave("Output_figures/Distribution_graph_with_median.pdf", width = 7, height = 7)


library(kableExtra)
descr_stat <- edgelist_indices_tidy %>% group_by(type) %>% 
  summarize(min = min(value), mean = mean(value), q3 = quantile(value, probs = 0.25),  max = max(value), 
            sd = sd(value), coef_var = sd/mean) 
descr_stat%>% 
  kbl( caption = "Distribution statistics of the indexes ",
       digits = 2, 
       col.names = c("Index","Min", "Mean", "Q3", "Max", "Standard Deviation", "Coefficient of Variation (Sd/Mean)")
       ) %>% 
  kable_styling(position ="center", 
                latex_options = c("HOLD_position"),
                font_size = 10) %>%
  # add_header_above(c("", "2018-2020" = 3, "2014-2020" = 1),
  #                  bold = T) %>%
  footnote(general = c("Source : Computed from the default probability matrix of InterVA4 (probbase Malaria : VL, HIV : VL)"
                       ), 
           general_title = "", threeparttable = T)%>%
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  kable_styling(latex_options = "hold_position") #%>%
# save_kable("Output_figures/Descriptive_statistics_index.pdf")
  





## b. Correlation ----
#install.packages('GGally')
library(ggplot2)
library(GGally)

ggpairs(edgelist_indices, columns = c(3,7:9),
        columnLabels = c("Absolut Norm", "Euclidean norm", "Normalised scalar product",
                         "Uniform normalisation"))
# ggsave("Output_figures/Correlogram.pdf", width = 7, height = 7)



## c. Heatmaps -----

# Preparing the data
M_indice_euclidien <- as.matrix(M_indice_euclidien) #Change into numerical values
M_indice_norme1 <- as.matrix(M_indice_norme1)
M_indice_scalar_product <- as.matrix(M_indice_scalar_product)

# For heatmaps more easily readable, we transform the 0 (from the top half diagonal of the matrix
# that are not calculated as they are symmetrical to the bottom ones into 1, to not visualise that part)
M_indice_euclidien[M_indice_euclidien == 0] <- 1
M_indice_norme1[M_indice_norme1== 0] <- 1
M_indice_scalar_product[M_indice_scalar_product == 0] <- 1



## setting a color palette
palf <- colorRampPalette(c("red", "white"), bias = 0.2) 

# Euclidean
heatmap(M_indice_euclidien, Rowv = NA, Colv = NA, col = palf(100),
        scale="none", margins=c(10,10),
        labCol = causetext$CAUSETXT.C.40, 
        labRow = causetext$CAUSETXT.C.40,
        main = "Heatmap of the Euclidean index of similarity")


# Norm 1 
heatmap(M_indice_norme1, Rowv = NA, Colv = NA, col = palf(100),
        scale="none", margins=c(10,10),
        labCol = causetext$CAUSETXT.C.40, 
        labRow = causetext$CAUSETXT.C.40,
        main = "Heatmap of the absolute norm index")

## Interactive heatmaps with plotly
library(plotly)

plot_ly(z = M_indice_euclidien, type = "heatmap", 
        x = causetext$CAUSETXT.C.40, 
        y = causetext$CAUSETXT.C.40, colors = palf(500),
        title="Heatmap of the Euclidean index of similarity")

plot_ly(z = M_indice_norme1, type = "heatmap",
        title="Heatmap of the absolute norm index of similarity",
        x = causetext$CAUSETXT.C.40, 
        y = causetext$CAUSETXT.C.40, colors = palf(500))

plot_ly(z = M_indice_scalar_product, type = "heatmap",
        title="Heatmap of the absolute norm index of similarity",
        x = causetext$CAUSETXT.C.40, 
        y = causetext$CAUSETXT.C.40, colors = palf(500))

plot_ly(z = M_indice_unif_normalisation, type = "heatmap",
        title="Heatmap of the absolute norm index of similarity",
        x = causetext$CAUSETXT.C.40, 
        y = causetext$CAUSETXT.C.40, colors = palf(500))



#### 4. Characteristic specific indexes ----


### a. HIV malaria prevalence adjusted indexes ----

### Probbase HH ----
## Download the probbase HH high HIV and high malaria prevalence

probbaseHH <- read.csv2("Input_probbase/HIV_malaria_prevalence_variations/probbaseHH_chiffres.csv",
                        row.names = 1, 
                        sep = ",")
probbaseHH$INDIC <- NULL #We don't need the indicator column, we eliminate it

## We redownload the list of indexes
edgelist_indexes <- read.csv("Output_indexes/VLVL_Edgelist_indexes.csv", row.name = 1)


## Create euclidean index with HH probbase
edgelist_euclidienHH <- as.data.frame(adjmat_to_edgelist(
  matricer(probbaseHH, indice_euclidien),
  undirected =FALSE,
  keep.isolates = TRUE
))

## We merge it with the general edgelist of indexes

edgelist_indexes_HH <- left_join(edgelist_indexes, edgelist_euclidienHH,
                              by = c("ego", "alter") ) %>%
  rename("indice_euclHH" = "value") 

# write.csv(edgelist_indexes_HH, "Output_indexes/HIV_malaria_variations/HH_VLVL_Matrix_Eucl_index_comparison.csv")


## Correlations ----

cor(edgelist_indexes_HH$indice_euclHH, edgelist_indexes_HH$indice_euclidien)

# We select more precisely associations with the two causes of death were their probabilities changed:
# HIV and malaria
edgelist_indexes_HH_HIV_malaria_CoD <- edgelist_indexes_HH %>% 
  filter(ego_name %in% c("HIV/AIDS related death", "Malaria") | 
         alter_name %in% c("HIV/AIDS related death", "Malaria"))

cor(edgelist_indexes_HH_HIV_malaria_CoD$indice_euclHH, 
         edgelist_indexes_HH_HIV_malaria_CoD$indice_euclidien)

## cor.test is for samples

### Probbase HL ----
## The only other possible difference is 

probbaseHL <- read.csv2("Input_probbase/HIV_malaria_prevalence_variations/probbaseHL_chiffres.csv",
                        row.names = 1, 
                        sep = ",")
probbaseHL$INDIC <- NULL #We don't need the indicator column, we eliminate it

## Create euclidean index with HH probbase
edgelist_euclidienHL <- as.data.frame(adjmat_to_edgelist(
  matricer(probbaseHL, indice_euclidien),
  undirected =FALSE,
  keep.isolates = TRUE
))

## See the difference for malaria HIV association
edgelist_Malaria_HIV <- left_join(edgelist_indexes, edgelist_euclidienHH,
                                 by = c("ego", "alter") ) %>%
  rename("indice_euclHL" = "value")%>% 
  filter(ego_name %in% c("HIV/AIDS related death", "Malaria"), 
           alter_name %in% c("HIV/AIDS related death", "Malaria"))


### b. Age-sex specific indexes ----

## To create an age-sex specific index we use the commented probbase specifying which indicators are not asked to certain sex or age groups.

probbase_comments <- read.csv("Input_probbase/probbaseVLVL_2012_comments.csv", sep = ";", header =T, row.names = 1) 
# We extract the column specifying the age sex groups not asked
probbase_comments_dontask <- probbase_comments %>% 
  select(INDIC.C.10, which(str_starts(names(probbase_comments), "DONTASK")))

probbase <- read.csv2("Input_probbase/probbaseVLVL_numbers.csv", sep = ",", row.names = 1)

edgelist_indexes <- read.csv("Output_indexes/VLVL_Edgelist_indexes.csv", row.name = 1)

## sex and age groups 

sex_names <- c("female", "male")
age_grp_names <- c("elder","midage","adult","child","under5","infant","neonate")


## function to create an age-sex specific matrix

#Test
# indic_name <- "INDIC.C.10"
# sex <- "male"
# age_grp <- "adult"

matricer_age_sex_indic <- function(age_grp,
                                   sex,
                                   probbase, ## M probbase matrix with indicator column
                                   indic_name, #name of the indicator column
                                   f ##  f the index function of choice
){  
  
  # Create an age and sex specific probbase matrix M
  age_sex_names <-c("elder","midage","adult","child","under5","infant","neonate","male","female")
  age_sex_to_discard <- age_sex_names[!(age_sex_names %in% c(age_grp, sex))]
  all_to_discard <- list(age_sex_to_discard)
  
  #Sekect all the indicators to discard
  for(i in c(2:ncol(probbase_comments_dontask))){
    
    indicators_to_discard <- list(probbase_comments_dontask$INDIC.C.10[probbase_comments_dontask[,i] %in% c(age_grp)|
                                                                         probbase_comments_dontask[,i] %in% c(sex)])
    
    all_to_discard[i] <- indicators_to_discard
    
    
  }
  all_to_discard_unlisted <- unlist(all_to_discard)
  
  M_indic <- probbase[!(probbase[[indic_name]] %in% all_to_discard_unlisted),] #matrix with indicator's names
  M <- M_indic[, -which(colnames(M_indic) == indic_name)] #numeric matrix
  
  # Extract only the age_sex_probabilities
  M_age_sex <- M[M_indic[[indic_name]] %in% c(age_grp, sex),]%>% mutate_all(as.numeric) %>% as.matrix()  #selecting age and sex probabilities, 
  # taking col number for age and sex from the matrix with indicators
  
  #calculate number of combinations of causes possible
  
  dimM<- as.matrix(dim(M))[2,1] ### store the number of columns i.e. possible causes
  
  comb <- combn(c(1:dimM), 2)# we compute all the possible combinations of columns
  dimcomb <- as.matrix(dim(comb))[2,1] # we put it in the right format 
  
  
  # We create an empty index matrix of the right size
  I <- matrix(data = 0, nrow = dimM, ncol = dimM)
  # We fill the matrix
  for (i in (1:dimcomb)){
    A <- comb[1,i] #Number associated with cause A
    B <- comb[2,i]
    
    #First we verify if there is probabilities = 0 in the age and sex categories
    
    A_B_age_sex_probabilities <- M_age_sex[, c(A,B)]
    presence_of_0 <- is.element(0,unlist(A_B_age_sex_probabilities))
    
    # We compute the indicator
    I[A,B] = case_when(presence_of_0 == T ~ 1,
                       T ~ f(unlist(M[,A]),unlist(M[,B])) 
    )
  }
  
  # rename the colums 
  colnames(I) <- colnames(M)
  rownames(I) <- colnames(M)
  
  return(I) # return our index matrix !!
}



## Calculate the number of indicators per age and sex combo 


count_age_sex_indic <- function(age_grp,
                                   sex,
                                   probbase, ## M probbase matrix with indicator column
                                   indic_name, #name of the indicator column
                                   f ##  f the index function of choice
){  
  
  # Create an age and sex specific probbase matrix M
  age_sex_names <-c("elder","midage","adult","child","under5","infant","neonate","male","female")
  age_sex_to_discard <- age_sex_names[!(age_sex_names %in% c(age_grp, sex))]
  all_to_discard <- list(age_sex_to_discard)
  
  #Sekect all the indicators to discard
  for(i in c(2:ncol(probbase_comments_dontask))){
    
    indicators_to_discard <- list(probbase_comments_dontask$INDIC.C.10[probbase_comments_dontask[,i] %in% c(age_grp)|
                                                                         probbase_comments_dontask[,i] %in% c(sex)])
    
    all_to_discard[i] <- indicators_to_discard
    
    
  }
  all_to_discard_unlisted <- unlist(all_to_discard)
  
  M_indic <- probbase[!(probbase[[indic_name]] %in% all_to_discard_unlisted),] #matrix with indicator's names
  
  return(nrow(M_indic))
}



## Compute per age and sex the number of indicators ----


nb_indicator_table <- data.frame(matrix(ncol = 3, nrow = length(sex_names)*length(age_grp_names)))
names(nb_indicator_table) <- c("sex", "age_grp", "nb_indic")

k = 0
for (i in (1:length(sex_names))){
  #Select sex
  sex <- sex_names[i]
  
  for (j in (1:length(age_grp_names))){
    # Determine row number
    k= k+1
    #Select age group
    age_grp <- age_grp_names[j]
    
    nb_indicator_table[k,]$sex <- sex
    nb_indicator_table[k,]$age_grp <- age_grp
    nb_indicator_table[k,]$nb_indic <- count_age_sex_indic(age_grp, sex, probbase,
                        "INDIC.C.10",
                        indice_euclidien)
    
  }
  
}

# write.csv(nb_indicator_table, "Output_figures/Age_sex_specific_index/Nb_indicator_age_sex.csv")


## Compute the similarity matrix for every combination ----



for (i in (1:length(sex_names))){
  sex <- sex_names[i]
  
  for (j in (1:length(age_grp_names))){
    
    age_grp <- age_grp_names[j]
    
    #print the age-sex combination
   print(sex)
   print(age_grp)
   
   #compute matrix
    M_age_sex <- matricer_age_sex_indic(age_grp, sex,
                           probbase, ## M probbase matrix with indicator column
                           "INDIC.C.10", #name of the indicator column
                           indice_euclidien ##  f the index function of choice
    )
    
  
    
    
    #save the matrix
    # write.csv(M_age_sex, paste0("Output_indexes/Age_sex_specific/VLVL_Matrix_Eucl_index_",sex,"_", age_grp, ".csv"))
    
    # plot matrix and print
    
    M_age_sex[M_age_sex == 0] <- 1
    
    png(paste0("Output_figures/Age_sex_specific_index/heatmap_eucl_index_",sex,"_", age_grp, ".png"),
        width = 8.27, height = 13,
        units = "in", res = 300)
    
    ## Setting plot margins
    print(heatmap(M_age_sex, Rowv = NA, Colv = NA, col = colorRampPalette(c("red", "white"), bias = 0.2)(100) ,
            scale="none", margins=c(10,10),
            labCol = causetext$CAUSETXT.C.40, 
            labRow = causetext$CAUSETXT.C.40,
            main = paste0("Heatmap of the Euclidean index of similarity - ", sex, " ", age_grp)
    ))                                         
    
    dev.off()
    
  }
    
  }
  

## Create an edgelist for comparison

for (i in (1:length(sex_names))){
  sex <- sex_names[i]
  
  for (j in (1:length(age_grp_names))){
    
    age_grp <- age_grp_names[j]
    
    #print the age-sex combination
    print(sex)
    print(age_grp)
    
    
    edgelist_index <- adjmat_to_edgelist( # the very useful function from netdiffuseR, 
      # allowing to transform an adjacency matrix to an edgelist
      matricer_age_sex_indic(age_grp,
                       sex,
                       probbase, ## M probbase matrix with indicator column
                       "INDIC.C.10", #name of the indicator column
                       indice_euclidien ##  f the index function of choice
      ),
      undirected =FALSE,
      keep.isolates = TRUE
    ) %>% as.data.frame() %>% 
      renommer(causetext, "ego", causetext$numero_cause
               , causetext$CAUSETXT.C.40, "ego_name") %>% 
      renommer(causetext, "alter", causetext$numero_cause
               , causetext$CAUSETXT.C.40, "alter_name") %>% 
      mutate(age_grp = age_grp, 
             sex = sex) %>% 
      filter(!(value == 1)) ## We take out the combinations that are not possible
    
    names(edgelist_index)[names(edgelist_index)=="value"] <- paste0("index_age_sex")
    
    #Put into dataframe
    edgelist_index <- left_join(edgelist_index, edgelist_indexes, ## We merge with the other indexes 
                                  by = c("ego", "alter", "time", "ego_name", "alter_name") )
    
    # we bind all the indexes together to create a comparison dataset
    if (i == 1 & j == 1) edge_list_index_comparison <- edgelist_index else edge_list_index_comparison <- bind_rows(edge_list_index_comparison, edgelist_index)
  
    }
  
}


# write.csv(edge_list_index_comparison, "Output_indexes/Age_sex_specific/EdgelistVLVL_age_sex_indexes_comparison.csv")

## Compute correlations ----


cor_table <- edge_list_index_comparison %>% group_by(sex, age_grp) %>% summarise(cor = round(cor(indice_euclidien, index_age_sex),3))
cor_table_nb_indicators <- left_join(nb_indicator_table, cor_table, by = c("sex", "age_grp"))

# write.csv(cor_table_nb_indicators, "Output_figures/Age_sex_specific_index/Cor_table.csv")


## Correlation graph ----

ggplot(edge_list_index_comparison, aes(x= indice_euclidien, y = index_age_sex)) +
  geom_point()+
  geom_smooth()+
  labs(x = reference_index,
       y = age_sex_sp_index)+
  facet_grid(rows = vars(age_grp),
    cols = vars(sex))


