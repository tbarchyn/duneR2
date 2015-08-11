# duneR - Thomas E. Barchyn
# Developed at the University of Calgary
# This software has no warranty whatsoever.

# Copyright Thomas E. Barchyn, 2015

# This file is part of duneR.

# duneR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# duneR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with duneR.  If not, see <http://www.gnu.org/licenses/>.


# PURPOSE : deal with and fix problems with outputs from degrib. Degrib
# is one of the better 'slicers' of grib files, but creates outputs that
# need some tidying up. The purpose of these functions is basically
# housecleaning.

assemble_degrib_collection <- function (collection) {
    # Function to assemble a collection of degrib outputs and do
    # necessary housecleaning. Returns a xts timeseries data frame
    # with attributes for variable names and units.
    
    # collection = a vector of names of degribbed files to be assembled.
    
    for (files in collection) {
        if (files == collection[1]) {
            # if this is the first file, start the degrib_collection object
            degrib_collection <- read_degrib (files)
        } else {
            # if this is subsequent files, read in the file, then merge
            degrib_collection_new <- read_degrib (files)
            
            # get the parameters and units from the files to be merged
            existing_parameters <- attributes(degrib_collection)$parameters
            existing_units <- attributes(degrib_collection)$units
            existing_parameters_new <- attributes(degrib_collection_new)$parameters
            existing_units_new <- attributes(degrib_collection_new)$units
            
            # merge
            degrib_collection <- merge (degrib_collection, degrib_collection_new)
            
            # fix the attributes
            xtsAttributes(degrib_collection) <- list (parameters = append(existing_parameters, existing_parameters_new))
            xtsAttributes(degrib_collection) <- list (units = append(existing_units, existing_units_new))
        }
        print (paste ('Completed read of file: ', files))
    }
    return (degrib_collection)
}

read_degrib <- function (filename) {
    # Function to read one degrib file and return an xts timeseries
    # with each of the constituent variables correctly referenced.
    
    # filename = the filename of the degribbed file
    
    # read in the table, set colnames and clean up times
    raw_degrib <- read.table (filename, skip = 1, sep = ',')
    names (raw_degrib) <- c('parameter', 'units', 'reftime', 'validtime', 'val')
    raw_degrib$units <- as.character (raw_degrib$units)
    raw_degrib$reftime <- format_degrib_time (raw_degrib$reftime)
    raw_degrib$validtime <- format_degrib_time (raw_degrib$validtime)

    # get the constituent levels
    raw_degrib_levels <- levels (raw_degrib$parameter)
    
    # for each level, cut out the data frame
    for (i in raw_degrib_levels) {
        # subset out the level
        subset_level <- raw_degrib[raw_degrib$parameter == i, ]
        units_spec <- subset_level$units[1]      # pull the units

        if (i == raw_degrib_levels[1]) {
            # if this is the first level, create the xts_degrib object
            xts_degrib <- format_degrib_xts (subset_level, i)
            xtsAttributes(xts_degrib) <- list (parameters = i, units = units_spec)
        } else {
            # if this is a subsequent level, merge with the existing xts_degrib object
            xts_degrib_new <- format_degrib_xts (subset_level, i)
            existing_parameters <- attributes(xts_degrib)$parameters
            existing_units <- attributes(xts_degrib)$units
            xts_degrib <- merge (xts_degrib, xts_degrib_new, all = T)
            xtsAttributes(xts_degrib) <- list (parameters = append(existing_parameters, i))
            xtsAttributes(xts_degrib) <- list (units = append(existing_units, units_spec))
        }
    }
    return (xts_degrib)
}

format_degrib_time <- function (time_col) {
    # Function to format a column of time, read directly from the
    # read.table function. Returns a POSIXlt formatted time code.
    
    # time_col = a vector of time inputs, this seems to read in
    #            as a numeric column.
    
    time_col <- as.character (time_col)
    time_col <- as.POSIXlt (time_col, format = "%Y%m%d%H%M", tz = "UTC")
    
    return (time_col)
}

format_degrib_xts <- function (subset_raw_degrib, lev_label) {
    # Function to take a subset of the raw_degrib associated
    # with a specific parameter code, and convert it to an indexed
    # xts timeseries. We can then match it up with other records confidently
    # in preparation for additional analyses. Returns indexed xts
    # dataframe.
    
    # subset_raw_degrib = a raw_degrib file
    # lev_label = name of the level to rename the values
    
    # first check for agreement between reftime and validtime. Throw warning
    # if there is disagreement. We will use the 'validtime' column though
    # for the index, and toss reftime if there is not agreement.
    nrow_disagreement <- nrow (subset_raw_degrib [subset_raw_degrib$reftime != subset_raw_degrib$validtime, ])
    if (nrow_disagreement > 0) {
        print (paste ('WARNING: disagreement between reftime and validtime in',
                        nrow_disagreement, 'records.'))
        # rename the column for reftime and coerce to character
        subset_raw_degrib$reftime <- as.character (subset_raw_degrib$reftime)
        names(subset_raw_degrib)[names(subset_raw_degrib) == 'reftime'] <- paste ('reftime', lev_label, sep = '_')
    } else {
        # if no disagreement, remove the reftime column
        subset_raw_degrib <- subset_raw_degrib [, names(subset_raw_degrib) != 'reftime']
    }
    
    # remove parameter and units cols
    subset_raw_degrib <- subset_raw_degrib [, names(subset_raw_degrib) != 'parameter']
    subset_raw_degrib <- subset_raw_degrib [, names(subset_raw_degrib) != 'units']
    
    # cut out the posixlt column from the dataframe
    time_index <- subset_raw_degrib$validtime
    
    # remake the data frame so it has a correct label
    subset_raw_degrib <- data.frame (val = subset_raw_degrib$val)
    names(subset_raw_degrib)[names(subset_raw_degrib) == 'val'] <- lev_label
        
    # create the xts data frame correctly indexed with the time_col external
    xts_degrib <- xts (x = subset_raw_degrib, order.by = time_index)
    return (xts_degrib)
}    

add_numeric_date <- function (xts_obj) {
    # Function to add numeric columns for the year, month, day, and hour to a
    # xts object. This makes monthly analyses, etc. easier.
    
    time_index <- index (xts_obj)
    
    # add year, month, day, hour columns
    xts_obj$year <- as.numeric(strftime(time_index, format = '%Y'))
    xts_obj$month <- as.numeric(strftime(time_index, format = '%m'))
    xts_obj$day <- as.numeric(strftime(time_index, format = '%d'))
    xts_obj$hour <- as.numeric(strftime(time_index, format = '%H'))
    return (xts_obj)
}
