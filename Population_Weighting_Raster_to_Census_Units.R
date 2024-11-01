library("terra")  # For raster data
library("sf")     # For vector data
library("plyr")
library("tigris") # For downloading census shapefiles
library("doBy")
library("tidyverse")
library("tidycensus")
library("lwgeom")
sf_use_s2(FALSE)  # S2 is for computing distances, areas, etc. on a SPHERE (using
                  # geographic coordinates, i.e., lat/lon in decimal-degrees); no need for this
                  # extra computational processing time if using PROJECTED coordinates,
                  # since these are already mapped to a flat surface. Here, PRISM
                  # is indeed in geographic coordinates, but the scale of areas we are 
                  # interested in is very small, and hence the error introduced by 
                  # ignoring the Earth's curvature over these tiny areas is negligible and
                  # a reasonable trade off given the dramatic reduction in processing time. Moreover,
                  # the areas we calculate are not an integral part of the process
                  # and any error in that step would not materially impact the final output
options(tigris_use_cache = TRUE)

# Check package version numbers
#
if (packageVersion("terra") < "1.5.34"   | packageVersion("sf") < "1.0.7" | 
    packageVersion("plyr")  < "1.8.7"    | packageVersion("tigris") < "2.0.4" |
    packageVersion("doBy")  < "4.6.19"   | packageVersion("tidyverse") < "1.3.1" |
    packageVersion("tidycensus") < "1.5" | packageVersion("lwgeom") < "0.2.8") {
  cat("WARNING: one or more packages are outdated. Please update packages to prevent potential errors. \n") }

# %%%%%%%%%%%%%%%%%%%%%%% USER-DEFINED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%% #

census_api_key("")    # Obtain an API key here: https://www.census.gov/data/developers.html
# NOTE: You should NOT save your API key within your code

input_data_dir <- "In_Dir/"   # Full pathway of the directory where your input data will be stored.
                              # ** If you change this, be sure to keep the final forward slash **
                              # "Input data" is any data set that is *not* the final, analytical dataset.

zip_file_dir <- "In_Dir/Zip_Files/" # Set the directory for where the .zip PRISM files
                                    # will be saved. If changing this directory, 
                                    # be sure to keep the final forward slash

output_data_dir <- "Out_Dir/" # Enter the full pathway of the directory where your output data will be stored.

stateFIPS <- "11"     # D.C. - for the purpose of this tutorial, we will only assess DC, a small area

year <- 2020          # For the purpose of this tutorial, we are calculating only 2020 data

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OVERVIEW  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# Step 1. Bulk-download PRISM data from FTP
# Step 2. Create PRISM extraction points
# Step 3. Calculate land-area weighted PRISM values at the census block level
# Step 4. Calculate *population*-weighted PRISM values at the census block group through county levels 
#
# This code is designed to illustrate the calculation of population-weighted mean
# values in PRISM for a single month in Washington, D.C.
#
# The code can be adapted to multiple times and locations, but requires substantial
# computing time. It is recommended that for nationwide analyses one use distributed 
# processing on a computing cluster, divided between time and space (for example,
# processing one month at a time for individual states).
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%% STEP 1: DOWNLOAD PRISM DATA FROM FTP %%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# This script downloads raw PRISM data at the ~4km resolution. Approximately ~800m
# resolution PRISM data is available from the PRISM Climate Group for a fee.
# Most applications do not need the enhanced resolution, and one should note
# that the higher-resolution product does not resolve urban heat islands; the modeling
# methodology does not incorporate urbanization characteristics such as LST, albedo,
# or imperviousness of surfaces (which are proxy indicators of UHI).
#
options(timeout=1000) # Set the max timeout (in seconds) for downloading files

# Identify the PRISM variables that you want to download using the syntax from
# the PRISM FTP. For the purpose of this tutorial, we are just using "tmax" and "tmin".
# Other options available: c("ppt", "tdmean", "tmean", "tmin", "vpdmax", "vpdmin")
#
vars <- c("tmax", "tmin")

# The baseline URL for the PRISM FTP directory for daily data. Note that
# other time steps are available as well. Explore the FTP site to find the
# relevant URL and syntax for the data sets that you need.
#
URL <- "https://ftp.prism.oregonstate.edu/daily/"

# Identify the years of data you need. NOTE: data < 6 months old are provisional.
# Use "stable" data whenever possible, and note that the filename will be different
# for provisional data (you will have to modify the code below accordingly).
#
years_data <- c(2020) # If downloading multiple years, enter range here, e.g., 2010:2020

for (i in 1:length(vars)) { 
  
  var <- vars[i]
  
  cat("---------------------------------------------------------------------\n")
  cat("Beginning download of", var, "PRISM data: variable", i, "of", length(vars), "\n")
  
  for (j in 1:length(years_data)) {
    
    year_data <- years_data[j]
    
    cat(".....Processing", year_data, "data \n")
    
    # Identify all of the days in that particular year in the format YYYYMMDD
    # This step is needed to account for Leap Days
    #
    days <- format(seq(as.Date(paste0(year_data, "-01-01")),
                       as.Date(paste0(year_data, "-12-31")), by = "days"),
                   format="%Y%m%d")
    
    for (k in 1:length(days)) {
      
      day <- days[k]
      
      dl_link <- paste0(URL, var, "/", year_data, "/PRISM_", var, "_stable_4kmD2_", days[k], "_bil.zip")
      dl_file <- paste0(zip_file_dir, "PRISM_", var, "_stable_4kmD2_", day, "_bil.zip")
      
      # Check to see if the file already exists; download if not
      #
      if (file.exists(dl_file)) {
        
        cat("Zip file for day", day, "already downloaded. Proceeding to next day. \n")
        
      } else {
        
        dl <- try(download.file(dl_link, destfile = dl_file))
        
        # Determine if the download failed
        #
        if (class(dl) == "try-error") {
          
          cat("ERROR! File did not download successfully for", day, "\n")
          cat("..... Pausing for 10 seconds and re-trying. \n")
          
          k <- 10
          while (k > 0) {
            
            cat("..... Attempt", ((10 - k) + 1), "out of 10... \n")
            Sys.sleep(10) # Pause for 10 seconds
            dl <- try(download.file(dl_link, destfile = dl_file))
            
            if (class(dl) != "try-error") { 
              
              cat("..... :) success! Moving on to next file \n"); break 
              
            } else {
              
              cat("..... :( unsuccessful! \n")
              k <- k - 1 
            }
          }
          
          if (k == 0) { cat("..... Error unresolved; file for", day, "has still not been downloaded. \n"); break }
          
        } else { cat(":) file for", day, "downloaded successfully \n") }   
      }
      
      # Unzip the downloaded files
      # Recommended to delete the original zip files manually after downloading is complete
      #
      unzip(dl_file, exdir = input_data_dir)
    }
  }
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%% STEP 2: CREATE PRISM EXTRACTION POINTS  %%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# %%%%%%%%%%%%%%%%%% STEP 2A. READ IN THE PRISM RASTER DATA %%%%%%%%%%%%%%%%%% #
#
prism_files <- list.files(input_data_dir, pattern = paste0(".*", "_", year, ".*bil$"))

# Automated QC check -- ensure all files are present
#
expected_days <- format(seq(as.Date(paste0(year, "-01-01")),
                            as.Date(paste0(year, "-12-31")), by = "days"), format = "%Y%m%d")
file_days <- substr(sapply(strsplit(prism_files, "_4kmD2_"), "[", 2), 1, 8)

if (length(prism_files) != length(expected_days) * length(vars)) {
  cat("ERROR: wrong number of files for", year, "(check the pattern syntax above) \n") } else { 
    if (length(which( !(expected_days %in%  file_days) )) > 0) {
      cat("ERROR: right number of files, but missing days (check duplicates) \n")
    } else { cat(":) all files are present \n") } }

#
# %%%%%%%%%%% STEP 2B. CREATE A FISHNET GRID OF THE RASTER EXTENT %%%%%%%%%%%% #
#
# Here, we are making a shapefile that is a fishnet grid of the raster extent.
# It will essentially be a polygon of lines surrounding each PRISM cell.
#
# Reference/credit: https://gis.stackexchange.com/a/243585
#
prism_files <- paste0(input_data_dir, prism_files)
prism_raster <- rast(prism_files[1])
prism_extent <- ext(prism_raster)
xmin <- prism_extent[1]
xmax <- prism_extent[2]
ymin <- prism_extent[3]
ymax <- prism_extent[4]

prism_matrix <- matrix(c(xmin, ymax,
                         xmax, ymax,
                         xmax, ymin,
                         xmin, ymin,
                         xmin, ymax), byrow = TRUE, ncol = 2) %>%
  list() %>% 
  st_polygon() %>% 
  st_sfc(., crs = st_crs(prism_raster))

# Create fishnet of the PRISM matrix
#
prism_rows <- dim(prism_raster)[1]
prism_cols <- dim(prism_raster)[2]
prism_fishnet <- st_make_grid(prism_matrix, n = c(prism_cols, prism_rows), 
                              crs = st_crs(prism_raster), what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))

# Automated QC check -- confirm same coordinate reference system (CRS) between
#                       the fishnet and PRISM raster
if ( !(all.equal(st_crs(prism_raster), st_crs(prism_fishnet))) ) {
  cat("ERROR: CRS's do not match \n") } else { cat(":) CRS's match \n") }

#
# %%%%%%%%%%%%%% STEP 2C. DOWNLOAD CENSUS BLOCK SHAPEFILE %%%%%%%%%%%%%%%%%%%% #
#
# In this step, we will download a shapefile of census blocks from the US Census
# Bureau. Alternatively, you can use your own block shapefile if you already have it
# saved locally.
#
blocks <- tigris::blocks(state = stateFIPS, year = year)

# Match the CRS of the blocks shapefile to the fishnet and PRISM data
#
blocks <- st_transform(blocks, crs = st_crs(prism_fishnet))

if (!isTRUE(all.equal(st_crs(blocks), st_crs(prism_fishnet)))) {
  cat("ERROR: CRS's don't match \n")  } else { cat(":) CRS's match \n") }

#
# %%%%%%%%%% STEP 2D. CREATE UNION BETWEEN FISHNET AND CENSUS BLOCKS %%%%%%%%% #
# Reference/credit: https://stackoverflow.com/a/68713743
#
my_union <- function(a,b) {
  st_agr(a) = "constant"
  st_agr(b) = "constant"
  op1 <- st_difference(a, st_union(b))
  op2 <- st_difference(b, st_union(a))
  op3 <- st_intersection(b, a)
  union <- rbind.fill(op1, op2, op3)
  return(st_as_sf(union))
}

#
# Subset fishnet to stateFIPS block so processing time is cut down
#
bbox <- st_bbox(blocks)
fishnet <- st_crop(prism_fishnet, bbox)

blocks <- st_make_valid(blocks) 
fishnet <- st_make_valid(fishnet)

# Create the union between the fishnet and blocks layers
#
fishnetblock <- my_union(fishnet, blocks)
fishnetblock$UniqueID <- 1:dim(fishnetblock)[1]

# Automated QC -- Check to see if the union has introduced any geometry errors
#                 and fix as appropriate
#
check <- try(st_make_valid(fishnetblock), silent = TRUE)

if (class(check)[1] == "try-error") {
  
  cat("There is an issue with the sf object \n")
  cat("..... Attempting fix \n")
  
  geo_types <- unique(as.character(st_geometry_type(fishnetblock, by_geometry = TRUE)))
  
  cat("..... Geometry types in sf object 'fishnetblock':", geo_types, "\n")
  
  for (j in 1:length(geo_types)) {
    
    fishnetblock_subset <- fishnetblock[which(st_geometry_type(fishnetblock, by_geometry = TRUE) == geo_types[j]),]
    if (j == 1) { updated_fishnetblock <- fishnetblock_subset; next }
    updated_fishnetblock <- rbind(updated_fishnetblock, fishnetblock_subset)
    
  }
  
  check2 <- try(st_make_valid(updated_fishnetblock), silent = TRUE)
  
  if (class(check2)[1] == "try-error") {
    cat("..... ERROR NOT RESOLVED \n") } else {
      cat("..... :) issue has been fixed! \n")
      
      updated_fishnetblock <- updated_fishnetblock[order(updated_fishnetblock$UniqueID),]
      if ( !(all.equal(updated_fishnetblock$UniqueID, fishnetblock$UniqueID)) ) {
        cat("ERROR: Unique IDs do not match \n") } else {
          cat(":) unique ID's match. Reassigning 'updated_fishnetblock' to 'fishnetblock' \n")
          fishnetblock <- updated_fishnetblock    
        }
    }
}


# Automated QC check -- ensure that there is a variable identifying the census
#                       geographies. Note that these variable names may change
#                       depending on the year of the census data.
#
state_id_var <- names(fishnetblock)[grep("^STATEFP", names(fishnetblock), ignore.case = TRUE)]
if (length(state_id_var) != 1) {
  cat("ERROR: missing variable name or multiple matches \n") } else {
    cat(":) variable name present. It is called", state_id_var, "\n") }

# Identify the polygons of the fishnet that do not intersect with the census
# data; drop them.
#
before_dim <- dim(fishnetblock)[1]
fishnetblock <- fishnetblock[which( !(is.na(fishnetblock[[state_id_var]])) ),]
after_dim <- dim(fishnetblock)[1]

cat("Dropped", before_dim - after_dim, "polygons that do not intersect with census data \n")

# Some polygons formed in the union are incredibly small -- this adds unnecessary
# computation time without materially reducing error. Drop the small polygons.
# NOTE: Typically, when calculating areas of polygons, you would want to convert to
#       a projected CRS appropriate for your study domain. For the purpose of identifying
#       negligibly small areas to drop here, the error introduced by using geographic
#       coordinates for calculating area at this scale is negligible.
#
fishnetblock$Area_m2 <- as.numeric(st_area(fishnetblock))

fishnetblock <- fishnetblock[which(fishnetblock$Area_m2 > 10),]

#
# %%%%%%%%%%%%%%%%%%% STEP 2E. CONVERT POLYGON TO POINTS %%%%%%%%%%%%%%%%%%%%% #
#
# The final step is to create the extraction points. This is a point shapefile
# that will enable us to extract PRISM data from an entire stack of rasters rather
# than individually processing zonal statistics on each raster layer. 
#
# NOTE: This step throws a warning message related to using geographic coordinates 
#       rather than a projected CRS. This step is only placing a point inside the 
#       polygon to identify which PRISM grid cell we need to extract from; as all
#       of the input data are on the same CRS and the spatial scale of the polygons
#       is extremely small, this does not introduce substantive error.
#
extraction_pts <- st_point_on_surface(fishnetblock)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%% STEP 3: CALCULATE BLOCK-LEVEL MEAN PRISM VALUES %%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# %%%%%%%%%%%%%%%%%% STEP 3A. DOWNLOAD POPULATION COUNTS %%%%%%%%%%%%%%%%%%%%% #
#
# Obtain block-level population data. Note that the Census Bureau only releases
# population counts at the block level in the Decennial Census (2000, 2010, 2020, etc.).
# Note also that block identifiers (i.e., GEOID or FIPS code) can change over time.
# When merging population data with non-Decennial Census years, you may need to 
# account for these changing identifiers.
#
# You can download Census data in a variety of ways. Here, we are using the
# tidycensus package, which uses the Census Bureau's API. You must first obtain 
# an API key. Visit https://www.census.gov/data/developers.html and select the
# large "Request a Key" image on the left side of the screen to request an API key.
# Enter your API key at the top of the script

# Block populations must be queried for each individual county. Identify all of the
# counties in the state being processed and then loop through them
#
counties <- tigris::counties(state = stateFIPS, year = year, cb = TRUE)
counties <- unique(counties[[grep("^COUNTYFP", names(counties), ignore.case = TRUE)]])

# Block populations are only available for Decennial Censuses
#
decennial_year <- ifelse(year %in% 2000:2009, 2000,
                         ifelse(year %in% 2010:2019, 2010,
                                ifelse(year %in% 2020:2029, 2020, NA)))

for (i in 1:length(counties)) {
  
  blockpop_bycounty <- get_decennial(geography = "block",
                                     variables = "P1_001N", # Use load_variables("pl", year = 2020) to see available vars
                                     year = decennial_year,
                                     state = stateFIPS,
                                     county = counties[i],
                                     sumfile = "pl")
  if (i == 1) { blockpop <- blockpop_bycounty; next }
  blockpop <- rbind(blockpop, blockpop_bycounty)
}


# Automated QC check -- confirm data is total population
#
if (length(which(blockpop$variable != "P1_001N")) > 0) {
  cat("ERROR: Not all pop. vars. \n") } else { cat(":) all pop. vars. \n") }

pop_var_name <- paste0("Pop_", year)
names(blockpop)[grep("^value$", names(blockpop), ignore.case = TRUE)] <- pop_var_name
blockpop[,c("NAME", "variable")] <- NULL

# Get the GEOID for the block shapefile and for the block populations file, which may 
# change depending on the year of Census data
#
GEOID_shapefile <- names(blocks)[grep("^GEOID", names(blocks), ignore.case = TRUE)]
GEOID_popfile <- names(blockpop)[grep("^GEOID", names(blockpop), ignore.case = TRUE)]

# Merge the population data with the shapefile
#
blocks <- merge(blocks, blockpop, by.x = GEOID_shapefile, by.y = GEOID_popfile, all.x = TRUE)
num_missing <- length(which(is.na(blocks[[pop_var_name]])))

if (num_missing > 0) {
  cat("WARNING:", num_missing, "blocks with missing population \n") 
} else { cat(":) no blocks with missing population \n") }

#
# %%%%%%%%%%%% STEP 3B. CALCULATE LAND-AREA WEIGHTED AVG BY BLOCK %%%%%%%%%%%% #
#
#
# Set the unique geographic identifier (GEOID) from the extraction_pts sf object
#
GEOID <- names(extraction_pts)[grep("^GEOID", names(extraction_pts))]

# Get the total area by block to calculate the spatial weight value (typically 1.0)
#
eqn <- as.formula(paste0("Area_m2 ~ ", GEOID))

# Define a function to calculate sums such that if all values are NA then it returns
# NA rather than 0.
#
sumfun <- function(x) { return(ifelse(all(is.na(x)), NA, sum(x, na.rm = TRUE))) }

ptstotal <- summaryBy(eqn, data = as.data.frame(extraction_pts), FUN = sumfun)

extraction_pts <- merge(extraction_pts, ptstotal, by = GEOID, all.x = TRUE)
extraction_pts$SpatWt <- extraction_pts$Area_m2 / extraction_pts$Area_m2.sumfun

# Convert the extraction_pts to a SpatVector object
#
extraction_pts <- terra::vect(extraction_pts)

#
# %%%%%%%%%%%% STEP 3C. EXTRACT PRISM VALUES FOR EACH BLOCK POINT %%%%%%%%%%%% #
# 
# In this step, we will use the extraction points to extract the PRISM grid cell
# underlying each portion of a census block across the entire raster stack of values.
#
# The PRISM data are organized as separate rasters by variable, so we need to loop
# through each of the raster stacks. 
#

for (i in 1:length(vars)) {
  
  cat("-------------------------------------------------------------------- \n")
  cat("Processing", vars[i], "\n")
  
  newvarname <- paste0(vars[i], "_C") # For naming the applicable column in the output data
  
  # Stack all of the daily files by year
  #
  prism_stack <- rast(prism_files[grep(vars[i], prism_files)])
  
  # Automated QC check -- confirm all days are present
  #
  if (dim(prism_stack)[3] != length(expected_days)) { 
    cat("ERROR: missing days in raster stack \n") } else { cat(":) all days present in raster stack \n") }
  
  prismpts <- terra::extract(prism_stack, extraction_pts)
  prismpts <- as_tibble(prismpts)
  prismpts <- bind_cols(as_tibble(extraction_pts[,c(GEOID, "SpatWt")]), prismpts) # The extraction is in the same order as the input extraction_pts
  prismpts$geometry <- NULL
  dates <- substr(sapply(strsplit(names(prismpts)[grep("PRISM_", names(prismpts))], "_4kmD2_"), "[", 2), 1, 8)
  names(prismpts)[grep("PRISM_", names(prismpts))] <- dates
  
  prismpts <- pivot_longer(prismpts, cols = all_of(dates), names_to = "PRISM_Date")
  names(prismpts)[which(names(prismpts) == "value")] <- newvarname
  
  # Before we calculate the final weighted average of Tmax and Tmin, we need to check for missing data.
  # If a value is NA on one of the polygons, then it will give an underestimate of the
  # temperature since the weights will no longer add to 1. Example: there are two
  # polygons, each with 50% area. If Tmax is 30 C in one and NA in the other, then
  # the area weighted average (which removes NA values) would give: (30 * 0.5) + (NA * 0.5) = 15 C.
  # Therefore, we need to re-weight the weights based on the availability of data.
  #
  eqn <- as.formula(paste0("SpatWt ~ ", GEOID, " + PRISM_Date"))
  avail <- summaryBy(eqn,
                     data = prismpts[which( !(is.na(prismpts[[newvarname]])) ),],
                     FUN = sumfun)
  
  # Merge this value back into the prismpts
  #
  prismpts <- merge(prismpts, avail, by = c(GEOID, "PRISM_Date"), all.x = TRUE)
  
  # Re-weight the area weight by dividing by total available weight
  #
  prismpts$SpatWt <- prismpts$SpatWt / prismpts$SpatWt.sumfun
  
  # Automated QC -- Check that the weights of *available data* all add to 1
  #
  eqn <- as.formula(paste0("SpatWt ~ ", GEOID, " + PRISM_Date"))
  check <- summaryBy(eqn,
                     data = prismpts[which( !(is.na(prismpts[[newvarname]])) ),],
                     FUN = sumfun)
  
  if (length(which(round(check$SpatWt.sumfun, 4) != 1)) > 0) {
    cat("ERROR: weights do not sum to 1 \n"); break 
  } else {
    cat(":) weights sum to 1 \n")
    prismpts$SpatWt.sumfun <- NULL
  }
  
  # Multiply the variable of interest (here "newvarname") by the weighting value and then
  # sum up the resultant values within census blocks. This is an area-weighted average;
  # we will use the same algorithm with population to get the population-weighted average
  # at block groups, tracts, and counties in Step 4 of the code.
  #
  tempvar <- paste0(newvarname, "_Wt")
  prismpts[[tempvar]] <- prismpts[[newvarname]] * prismpts[["SpatWt"]]
  
  eqn <- as.formula(paste0(tempvar, " ~ ", GEOID, " + PRISM_Date"))
  final <- summaryBy(eqn, data = prismpts, FUN = sumfun)
  
  # Automated QC -- Confirm that the dimensions are as expected
  #
  if (length(unique(extraction_pts[[GEOID]][[GEOID]])) * length(unique(prismpts$PRISM_Date)) != dim(final)[1]) {
    cat("ERROR: incorrect dimensions of final data frame \n"); break
  } else { cat(":) dimensions of final data frame are as expected \n") }
  
  # Merge back in the block-level populations
  #
  final <- merge(final, as.data.frame(blocks[,c(GEOID, pop_var_name)]),
                 by = GEOID, all.x = TRUE)
  
  names(final)[grep(paste0("^", vars[i]), names(final))] <- newvarname
  
  final <- final[,c(GEOID, "PRISM_Date", pop_var_name, newvarname)]
  
  if (i == 1) {
    prism_block <- final
    cat(".....Successfully processed", vars[i], "\n")
  } else {
    if (all.equal(prism_block[[GEOID]], final[[GEOID]])) {
      prism_block[[newvarname]] <- final[[newvarname]]
      cat(".....Successfully processed", vars[i], "\n")
    } else {
      cat("ERROR: GEOID mismatch \n"); break
    }
  }
  
  if (i == 1) { final_output <- final; next }
  if (all.equal(final_output[[GEOID]], final[[GEOID]]) &
      all.equal(final_output[["PRISM_Date"]], final[["PRISM_Date"]])) {
    final_output <- cbind(final_output, final[[newvarname]])
    names(final_output)[length(final_output)] <- newvarname
  } else { cat("ERROR: GEOID and/or Date do not match between variables! \n"); break }
  
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%% STEP 4: AGGREGATE UP FROM BLOCK TO BG, TRACT, AND COUNTY %%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#
# The final step is to calculate population-weighted averages by variable at the
# levels of block groups (bg), census tracts, and counties. Census blocks nest
# within every other Census boundary; however, the only units in which they nest
# based on their GEOID (i.e., FIPS code) are block group, tract, county, and state.
# Other census units use different administrative boundaries. To aggregate up
# to these units, e.g., zip code tabulation area (ZCTA), you will need to
# use a crosswalk to link each census block to its accompanying ZCTA.
#
# By contrast, the GEOID variable contains the unique identifiers in which the
# block is nested:
#
# GEOID:  110010001011000 = block     (15 characters [+ 1 letter sometimes in intra-decadal years])
#         110010001011 = block group  (12 characters)
#         11001000101 = tract         (11 characters)
#         11001 = county              (5 characters)
#         11 = state                  (2 characters)
#
# Automated QC: check for variable type of GEOID
#               Common error: GEOID is read in as integer/numeric and the leading 0 gets dropped
#               for GEOIDs with state FIPS < 10
#
if (is(prism_block[[GEOID]])[1] != "character") { print("ERROR: wrong variable type for GEOID") 
} else { print(":) GEOID is correct type") }

# Create census ID variables for block group, tract, and county by 
# pulling out the first 12, 11, and 5 characters, respectively
#
prism_block$GEOID_BG <- substr(prism_block[[GEOID]], 1, 12)
prism_block$GEOID_Tract <- substr(prism_block[[GEOID]], 1, 11)
prism_block$GEOID_County <- substr(prism_block[[GEOID]], 1, 5)

# List relevant variable names
#
varnames <- paste0(vars, "_C") 

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%% CALCULATE POPULATION WEIGHTED MEANS BY CENSUS GEOGRAPHY %%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

censusgeos <- c("BG", "Tract", "County")
subgeoid <- c("GEOID_BG", "GEOID_Tract", "GEOID_County")

for (j in 1:length(censusgeos)) {
  
  cat("-----------------------------------------------------------\n")
  cat("Processing", censusgeos[j], "\n")
  
  for (i in 1:length(varnames)) {
    
    cat("..........Processing variable", varnames[i], "\n")
    
    # Calculate total population by geography (denominator of spatial weight)
    #
    pop_subgeo <- summaryBy(as.formula(paste0(pop_var_name, " ~ ", subgeoid[j], " + PRISM_Date")),
                            data = prism_block, FUN = sumfun)
    subgeopopvar <- names(pop_subgeo)[grep("^Pop", names(pop_subgeo))]
    
    # Merge these block group populations with the full block dataset
    #
    prism_block_subgeo <- merge(prism_block, pop_subgeo, by = c(subgeoid[j], "PRISM_Date"), all.x = TRUE)
    
    # Calculate weight for blocks within block groups
    #
    subgeowt <- paste0(censusgeos[j], "Wt")
    prism_block_subgeo[[subgeowt]] <- prism_block_subgeo[[pop_var_name]] / prism_block_subgeo[[subgeopopvar]]
    
    # Calculate available weights, which accounts for missingness in eventual output
    #
    avail <- summaryBy(as.formula(paste0(subgeowt, " ~ PRISM_Date + ", subgeoid[j])),
                       data = prism_block_subgeo[which( !(is.na(prism_block_subgeo[varnames[i]])) ),],
                       FUN = sumfun)
    
    # Merge available weight with raw weight
    #
    prism_block_subgeo <- merge(prism_block_subgeo, avail, by = c(subgeoid[j], "PRISM_Date"), all.x = TRUE)
    
    # Mark as NA any block-group-day/tract-day/county-day in which <50% of data are available
    #
    availwt <- paste0(subgeowt, ".sumfun")
    prism_block_subgeo[[availwt]][which(prism_block_subgeo[[availwt]] < 0.50)] <- NA
    
    # Adjust weight based on availability of data
    #
    prism_block_subgeo[[subgeowt]] <- prism_block_subgeo[[subgeowt]] / prism_block_subgeo[[paste0(subgeowt, ".sumfun")]]
    
    varwt <- paste0(varnames[i], "_Wt")
    prism_block_subgeo[[varwt]] <- prism_block_subgeo[[varnames[i]]] * prism_block_subgeo[[subgeowt]]
    
    # Automated QC -- Check to see if the weights all summed correctly
    #
    check1 <- summaryBy(as.formula(paste0(subgeowt, " ~ ", subgeoid[j], " + Date_PRISM")),
                        data = prism_block_subgeo[which( !(is.na(prism_block_subgeo[varnames[i]])) ),],
                        FUN = sumfun)
    
    numdays <- length(seq(min(as.Date(pop_subgeo$PRISM_Date, format = "%Y%m%d")),
                          max(as.Date(pop_subgeo$PRISM_Date, format = "%Y%m%d")),
                          by = "days"))
    if (length(which(round(check1[[paste0(subgeowt, ".sumfun")]], 4) != numdays)) > 0 |
        length(which(is.na(check1[[paste0(subgeowt, ".sumfun")]]))) > 0) {
      cat("WARNING: weights do not all sum to 1 \n") 
      cat("...... total number of rows marked as NA:", length(which(is.na(prism_block_subgeo[[varwt]]))), "\n")
      allNA <- check1[[subgeoid[j]]][which(is.na(check1[[paste0(subgeowt, ".sumfun")]]))]
      cat("...... number of geographies with NA on all days:", length(allNA), "\n")
      missingNon0 <- length(which(prism_block_subgeo[[pop_var_name]][which(prism_block_subgeo[[subgeoid[j]]] %in% allNA)] != 0))
      if (missingNon0 > 0) { cat("...... ERROR: some non-zero pop blocks are missing ALL days \n"); break }
      cat("............ >> among these,", missingNon0, "have non-zero populations \n")
      cat("...... among non-zero pop. geos, total number of geo-days missing:", sumfun(numdays - round(check1[[paste0(subgeowt, ".sumfun")]], 4)), "\n")
    } else { cat(":) all weights sum to 1 w/ no NA values \n") }
    
    # Sum the partial weights of the variable to get the final population-weighted value
    #
    final_wt <- summaryBy(as.formula(paste0(varwt, " ~ PRISM_Date + ", subgeoid[j])),
                          data = prism_block_subgeo, FUN = sumfun)
    
    names(final_wt)[grep(varnames[i], names(final_wt))] <- varnames[i]
    
    if (i == 1) { final <- final_wt; next }
    
    if (all.equal(final[[subgeoid[j]]], final_wt[[subgeoid[j]]])) {
      final <- cbind(final, final_wt[[varnames[i]]])
      names(final)[length(final)] <- varnames[i]
    } else { cat("ERROR: GEOID mismatch \n") }
    
  } # end of variable loop
  
  # Rename the "final" df to indicate which census geography it is
  #
  assign(paste0("final_", censusgeos[j]), final)
}

saveRDS(final_BG, paste0(output_data_dir, "PRISM_PopWt_", year, "_", stateFIPS, "_BlockGroup.Rds"))
saveRDS(final_Tract, paste0(output_data_dir, "PRISM_PopWt_", year, "_", stateFIPS, "_Tract.Rds"))
saveRDS(final_County, paste0(output_data_dir, "PRISM_PopWt_", year, "_", stateFIPS, "_County.Rds"))
