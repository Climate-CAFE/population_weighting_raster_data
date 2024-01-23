# Population-Weighting PRISM Data

This document describes the basics for population-weighting gridded meteorological data from the [PRISM Climate Group](https://www.prism.oregonstate.edu/).
Although this syntax uses PRISM data as an example, the methodology can be applied similarly to any gridded (i.e., raster) data. The script is completely self-contained, meaning that all the data needed to process the tutorial will be downloaded directly within the code.

## Why Use Population Weighting?

Imagine a large, hypothetical county that consists of vastly different terrain, with one half consisting of high-elevation mountains and the other half with low-elevation plains. Suppose also that the entire population lives in the low-elevation plains. If we are interested in population exposures to air temperature, then we want to know what the air temperature is in the county <i>where people live</i>. In this hypothetical example, the mountainous part of the county (which would have lower air temperatures on average) would have weights of 0 in the population-weighted average and hence not be factored in to the calculation. The result would be an average temperature that is more reflective of what people in the county actually experienced than we would get if we had taken a simple average across the entire county.

Note: Population weighting should only be used when the county-level estimates are being used in the context of human populations. 

## Overview of Process

There are many ways to calculate population-weighted averages. The method described here balances computational complexity with data validity (see Limitations). There are three overarching steps:
1. Create PRISM extraction points
    + This step conducts a union between a shapefile of census blocks and a fishnet of the PRISM grid and then converts it to points
    + The idea is that each point represents a portion of a census block (the smallest administrative unit) that falls within a respective PRISM grid cell
3. Calculate land-area weighted PRISM values at the census block level
    + This step uses area weighting because the census block is already the smallest areal unit for which the Census Bureau provides population counts
    + If a block falls 50% in one PRISM grid cell and 50% in another grid cell, then the value assigned to that block is a simple average of the two
5. Calculate *population-weighted* PRISM values for larger administrative units
    + Census blocks nest within every other administative unit provided by census
    + By merging the block-level PRISM values with the block-level populations, we can then aggregate up to any other unit by calculating an average weighted by the block population

 ## Limitations

 + This tutorial includes population weighting from census blocks to block groups, tracts, and counties, which can be easily achieved since these units nest neatly within one another and follow a predictable GEOID naming convention: block GEOIDs are 15-character identifiers where the first 12 characters represent the block group, the first 11 characters represent the census tract, and the first 5 characters represent the county. Although blocks nest within ZCTAs, they have separate GEOIDs; to aggregate up from block to ZCTA, a separate process must be followed that involves a relationship file between blocks and ZCTAs.
 + The spatial size of census blocks depends on population density, so blocks in rural areas are much larger than those in urban areas. Since step 2 of the overall process is an area-weighted average, this means that census areas (e.g., ZCTAs, tracts, etc.) that contain large, sparsely populated blocks may have more exposure-measurement error in their population-weighted averages than urban areas, where blocks frequently fall completely within a single PRISM grid cell. Higher-resolution population estimates exist that can be used for more highly resolved population weighting, if needed for your application.
 + Census only provides population counts at the block level in the Decennial Census products (2000, 2010, 2020, etc.). Although estimated adjustments to the block-level populations can be made using ACS data and some assumptions about the stationarity of within-block-group spatial distributions of population, those methods have not been implemented here. Therefore, when using the Decennial Census counts to calculate population-weighted exposures over the course of the decade, the exposure-measurement error of the population-weighted averages may increase over time as populations change within the decade.
 + Census boundaries change over time: substantial changes occur to blocks, block groups, and census tracts at every Decennial Census (2000, 2010, 2020, etc.), and additional changes can occur between the decades as well. This includes changes to county boundaries or names, as well as small corrections and adjustments made to administrative boundaries. Ensure that you are using an appropriate block-level shapefile for the time period of your data, and keep in mind that geographic identifiers (i.e., GEOID or FIPS codes) are subject to change over time.
