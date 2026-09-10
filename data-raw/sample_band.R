## code to prepare `sample_band` dataset goes here

bandpath <- "C:/Users/kdybala/Documents/R_projects/palodatavis/rawdat/allnumb.dbf"  # local copy of allnumb
band_raw = foreign::read.dbf(bandpath) |>  
  dplyr::select(.data$INITIALS:.data$COM) 
  # slow because this is a large database!
str(band_raw)

# take a sample of a few years worth, with multiple locations
sample_band = band_raw |> 
  filter(DATE >= '2009-01-01' & DATE <= '2012-12-31')

usethis::use_data(sample_band, overwrite = TRUE)
