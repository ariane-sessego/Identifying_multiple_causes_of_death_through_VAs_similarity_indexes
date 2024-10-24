##### Computing the individual-based index of similarity based on InterVA-4 ####

#### Script for the replicability of the resluts presented in
### "Studying multiple causes of death in the absence of death certificates:
### taking advantage of verbal autopsies"

#### Ariane SESSEGO under the supervision of Geraldine DUTHE (INED, France), 
#### in collaboration with Bruno Lankoandé and Dianou Kassoum (ISSP, Ouagadougou) 
#### March 2022 

#### R.version : 4.1.2

## Preparation -----
## Packages -----
library(questionr)
library(tidyverse)
library(FactoMineR)
library(haven)
library(dplyr)
library(xlsx)
library(igraph)
library(ggplot2)
library(kableExtra)
library(rJava)
library(stringr)
library(lubridate)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

## Importing the example database - a simulated example database, as the data cannot be shared.
detailed_VA <- read.csv2("Input_VA_example/Detailed_VA_examples.csv",
                      sep = ",")
detailed_VA$X <- NULL

## Probbase - importing the probbase with correct level of HIV and malaria prevalence
probbase <- read.csv("Input_probbase/probbaseVLVL_numbers.csv")
probbase$X <- NULL

## Causetext
library(xlsx)
causetext <-  read.csv("Input_probbase/causetext_AS.csv", sep = ",")
causetext <- causetext[4:63,]

## 1. Creating the function -----

## a. Defining the functions of the indexes  -----

# Let A and B be vectors that associate the probability of
# presenting each symptom given cause a and b respectively

## The Euclidean index
indice_euclidien <- function(A,B){
  A <- as.matrix(as.numeric(A)) # Transforming the vectors into numeric vectors to allox for computation  
  B <- as.matrix(as.numeric(B))
  indice <- norm(A - B, type = "F")/ norm(A + B, type = "F") # type = "F" i.e. chosing the Euclidean norm
  return(indice)
}

# The absolute norm index
indice_norme1 <- function(A,B){
  A <- as.matrix(as.numeric(A)) # Transforming the vectors into numeric vectors to allox for computation  
  B <- as.matrix(as.numeric(B))
  indice <- norm(A - B, type = "O")/ norm(A + B, type = "O") # type = "0" i.e. chosing norm 1
  return(indice)
}

#### b. recoding causes into numbers -----

# To facilitate the computation we recode causes into their identifying number according to causetext

# function to rename
renommer <- function(
    base_var, # base o??? l'on veut renommer les modalit???s d'une variable
    base_noms # base avec les equivalences des noms 
    , arenommer_entreguill # base 
    , equivalent, nv_nom
    , nvnom_guill #nom de la nouvelle variable cr??????e
){
  equivalent <- as.character(equivalent)
  # on extrait les deux var qui nous int???resse puis on fait une fusion
  C <- bind_cols(equivalent, nv_nom)
  # on renome pour que ce soit plus simple
  colnames(C) <- bind_cols(arenommer_entreguill, nvnom_guill)
  # on renomme pour que les cl???s aient le m???me nom pq apparemment il veut pas le faire autrement
  
  base_var <- left_join(base_var, C, by = arenommer_entreguill )
  
  return(base_var)
}

# We create a variable cause1_n with the number corresponding to cause 1
detailed_VA$cause1 <- as.character(detailed_VA$cause1)
detailed_VA<- renommer(detailed_VA, causetext, "cause1", causetext$V3, causetext$numero_cause, "cause1_n")

# We do the same for cause 2
detailed_VA$cause2 <- as.character(detailed_VA$cause2)
detailed_VA<- renommer(detailed_VA, causetext, "cause2", causetext$V3, causetext$numero_cause, "cause2_n")



#### c. Function to compute the index given a specific row (= one VA)-----



indice_empirique <- function(individu, # row
                             probbase, #probbase used
                             fonction_indice){ #the index function of choice
  
  ## Isolating the positively declared symptoms
  symptomes <- select(individu, 2:246) # Selecting all the symptoms used by the InterVA algorithm
  symptomes <- as.numeric(symptomes) # into a vector
  symptomes <- c(1, symptomes) # Add a 1 at the beginning to allow to take into account the a priori prevalence of the cause
  
  ## Extracting the vactor of probability of cause 1
  causeA <- individu$cause1_n  #  we identify the 1 cause through its number (equivalence name number in causetext)
  probaA <- probbase[,causeA+1]  # we extract the probability vector of cause 1 
  # causeA+1 because column 1 is the name of the symptoms
  
  # We multiply the two vectors to have the probability vector given the reported symptoms
  A <- symptomes*probaA
  
  ## We do the same for cause 2
  causeB <- individu$cause2_n  # we identify the second cause through its number (equivalence name number in causetext)
  probaB <- probbase[,causeB+1] 
  
  B <- symptomes*probaB
  
  ### We compute the index
  indice <- fonction_indice(A = A,B = B)
  
  return(indice)
}


##Example computing the index on one VA
indice_empirique(detailed_VA[1,], probbase, indice_euclidien)


## We calculate the empirical index for the whole dataset :

for (i in nrow(detailed_VA)){
  #with euclidean norm
  detailed_VA$indice_empirique_euclidien[i] <- indice_empirique(detailed_VA[i,], probbase, indice_euclidien)
  # with absolute norm
  detailed_VA$indice_empirique_norme1[i] <- indice_empirique(detailed_VA[i,], probbase, indice_norme1)
  
}

