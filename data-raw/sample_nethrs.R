## code to prepare `nethrs_sample` dataset goes here

nethrspath <- "C:/Users/kdybala/Documents/R_projects/palodatavis/rawdat/allpalonthrs.dbf" # local copy of entire database
effort_raw <- foreign::read.dbf(nethrspath) |> select(PROJECT:DUPE) 
unique(effort_raw$LOCATION)
# PN, PGUP, MUHO, PIGU, LACR, RECR, G5, CT, PT, HUMP, PEXA, PALO
summary(effort_raw$DATE) #1976-01-03 through 2025-12-31

# take a sample of a few years worth, with multiple locations
sample_nethrs = effort_raw |> 
  filter(DATE >= '2009-01-01' & DATE <= '2012-12-31')

usethis::use_data(sample_nethrs, overwrite = TRUE)
