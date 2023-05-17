# Creating a survival-swimmer plot in R(http://rstudio-pubs-static.s3.amazonaws.com/389060_a0ea812d8f8d490393bd1c65a6dcddef.html)
# Importing the required libraries
rm(list = ls())
library(magrittr)
library(stringi)
library(readr)   # Reading in the dataset
library(ggplot2) # Viewing the dataset
library(forcats) # Sorting factors
library(RColorBrewer) # Plot colours
library(dplyr, warn.conflicts=FALSE)   # Manipulating the dataframes
library(purrr, warn.conflicts=FALSE)   # Manipulating dataframe metadata
library(zoo, warn.conflicts=FALSE)     # Filling in  NA values
library(reshape2) # Reformmating dataframes 

# Importing the dataset
swimmer_file = "https://blogs.sas.com/content/graphicallyspeaking/files/2014/06/Swimmer_93.txt"
col.names = c("subjectID", "stage", "startTime", "endTime", 
              "isContinued", "responseType", "responseStartTime", "responseEndTime", "Durable")
df <- readr::read_lines(swimmer_file) %>%
  # Split by line recursion (\r\n)
  stringi::stri_split(fixed="\r\n", simplify=TRUE) %>%
  # Take only lines starting with a number (sample id)
  .[grepl("^[0-9]+", .)] %>%
  # Remove spaces from response column
  gsub(pattern="\\sresponse", replacement="_response") %>%
  # Remove spaces from stage column
  gsub(pattern="Stage\\s",  replacement="Stage_") %>%
  # Some lines missing 'Stage' and 'isContinued' column. 
  # Replace any set of 8 or more spaces with ' . '
  gsub(pattern="\\s{8,}", replacement=' . ') %>%
  # Split strings by spaces, do not include empty strings as columns
  stringi::stri_split(fixed=" ", simplify=TRUE, omit_empty=TRUE) %>%
  # Convert to dataframe
  as.data.frame(stringsAsFactors=FALSE) %>%
  # Set the column names
  purrr::set_names(col.names) %>%
  # We need to do some more cleaning up of the dataframe
  # Convert all . to NAs
  dplyr::na_if(".") %>%
  # Fill NAs in Stage column
  dplyr::mutate(stage=zoo::na.locf(stage)) %>%
  # Turn isContinued into boolean
  dplyr::mutate(isContinued=dplyr::if_else(isContinued=="FilledArrow", TRUE, FALSE, missing=FALSE)) %>%
  # Convert stage variable to factor, remove underscore
  dplyr::mutate(stage = as.factor(gsub(pattern="_", replacement=" ", x=stage))) %>%
  # Remove underscore from response types 
  dplyr::mutate(responseType = gsub("_", " ", responseType)) %>%
  # Change Durable from character to numeric
  dplyr::mutate(Durable = as.numeric(Durable)) %>%
  # Change Time variables from character to numeric
  dplyr::mutate_at(vars(dplyr::ends_with("Time")), as.numeric)

# Viewing the data.
# Let’s have a look at the data using the glimpse tool. We can see that the data frame is ‘tidy’. There is one row for every observation. A tidy dataframe may mean there are multiple entries for a given subject.
df %>% dplyr::glimpse()













