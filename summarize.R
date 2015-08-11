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

# PURPOSE: this file contains functions for summarizing timeseries data
# in different manners (e.g., monthly, yearly, etc.)

summarize_yearly <- function (x, years, summary_type, flux_col, record_duration_col,
                                wind_azimuth_col, ...) {
    # Function to compute yearly summaries of data. The year argument
    # is used as a factor level. The type argument is used to denote the
    # type of summary required. Some of these summaries could be rather
    # dune specific, the rationale behind doing this transparently. The
    # function returns a dataframe for each year, with named columns.
    
    # Yearly summaries are a meaningful time series technique given the
    # pronounced seasonality in weather data.
    
    # x = a dataframe, vector, or matrix
    # years = a vector of year values
    # summary_type = the type of summary required
    # flux_col = set flux column name for resultant calculations
    # record_duration_col = set record_duration column name for resultant calculations
    # wind_azimuth_col = set wind azimuth column name for resultant calculations    
    # ... = any other arguments passed down to the individual summary_type
    
    # deal with optional arguments
    if (missing (flux_col)) {
        flux_col <- ''
    }
    if (missing (record_duration_col)) {
        record_duration_col <- ''
    }
    if (missing (wind_azimuth_col)) {
        wind_azimuth_col <- ''
    }
    
    years <- factor (years)                         # factorize the years
    oput <- data.frame (year = levels(years))       # output dataframe

    years_num <- as.numeric(as.character(years))    # make numeric years
    oput$year <- as.numeric(as.character(oput$year))# force numeric
    x <- data.frame (x)                             # ensure this is a dataframe

    # create a prototype dataframe that will be cbinded to oput
    prototype <- data.frame (val = rep(NA, length (levels(years))))

    # loop through the names
    for (n in names(x)) {
        # calculate means
        if (summary_type == 'mean') {
            res <- prototype
            for (i in 1:nrow(res)) {
                res$val[i] <- mean (x[years_num == oput$year[i], n], ...)
            }
            names(res) <- paste(n, '_mean', sep = '')
            oput <- cbind(oput, res)
        }
        
        # calculate quantiles
        if (summary_type == 'quantiles') {
            # set the probabilities to report
            probabilities <- c(0.1, 0.25, 0.5, 0.75, 0.9)
            
            # loop through the probabilities
            for (p in probabilities) {
                res <- prototype
                for (i in 1:nrow(res)) {
                    res$val[i] <- quantile (x[years_num == oput$year[i], n], probs = p, ...)
                }
                names(res) <- paste(n, '_q', p, sep = '')
                oput <- cbind(oput, res)
            }
        }
        
        # sum
        if (summary_type == 'sum') {
            res <- prototype
            for (i in 1:nrow(res)) {
                res$val[i] <- sum (x[years_num == oput$year[i], n], ...)
            }
            names(res) <- paste(n, '_sum', sep = '')
            oput <- cbind(oput, res)
        }
    }
    
    # compute resultant
    if (summary_type == 'resultant') {
        oput$az <- NA
        oput$mag <- NA
        
        for (i in 1:nrow(oput)) {
            x_sub <- x[years_num == oput$year[i], ]
            res <- calc_resultant (flux = x_sub[, flux_col],
                                          record_duration = x_sub[, record_duration_col],
                                          wind_azimuth = x_sub[, wind_azimuth_col] )
            oput$az[i] <- res$az
            oput$mag[i] <- res$mag
        }
    }
    
    return (oput)
}  

summarize_monthly <- function (x, month, summary_type, flux_col, record_duration_col,
                                wind_azimuth_col, ...) {
    # Function to compute monthly summaries of data. The month argument
    # is used as a factor level. The type argument is used to denote the
    # type of summary required. Some of these summaries could be rather
    # dune specific, the rationale behind doing this transparently. The
    # function returns a dataframe for each month, with named columns.
    
    # monthly summaries are a meaningful time series technique given the
    # pronounced seasonality in weather data.
    
    # x = a dataframe, vector, or matrix
    # month = a vector of month values
    # summary_type = the type of summary required
    # flux_col = set flux column name for resultant calculations
    # record_duration_col = set record_duration column name for resultant calculations
    # wind_azimuth_col = set wind azimuth column name for resultant calculations    
    # ... = any other arguments passed down to the individual summary_type
    
    # deal with optional arguments
    if (missing (flux_col)) {
        flux_col <- ''
    }
    if (missing (record_duration_col)) {
        record_duration_col <- ''
    }
    if (missing (wind_azimuth_col)) {
        wind_azimuth_col <- ''
    }
    
    month <- factor (month)                         # factorize the month
    oput <- data.frame (month = levels(month))       # output dataframe

    month_num <- as.numeric(as.character(month))    # make numeric month
    oput$month <- as.numeric(as.character(oput$month))# force numeric
    x <- data.frame (x)                             # ensure this is a dataframe

    # create a prototype dataframe that will be cbinded to oput
    prototype <- data.frame (val = rep(NA, length (levels(month))))

    # loop through the names
    for (n in names(x)) {
        # calculate means
        if (summary_type == 'mean') {
            res <- prototype
            for (i in 1:nrow(res)) {
                res$val[i] <- mean (x[month_num == oput$month[i], n], ...)
            }
            names(res) <- paste(n, '_mean', sep = '')
            oput <- cbind(oput, res)
        }
        
        # calculate quantiles
        if (summary_type == 'quantiles') {
            # set the probabilities to report
            probabilities <- c(0.1, 0.25, 0.5, 0.75, 0.9)
            
            # loop through the probabilities
            for (p in probabilities) {
                res <- prototype
                for (i in 1:nrow(res)) {
                    res$val[i] <- quantile (x[month_num == oput$month[i], n], probs = p, ...)
                }
                names(res) <- paste(n, '_q', p, sep = '')
                oput <- cbind(oput, res)
            }
        }
        
        # sum
        if (summary_type == 'sum') {
            res <- prototype
            for (i in 1:nrow(res)) {
                res$val[i] <- sum (x[month_num == oput$month[i], n], ...)
            }
            names(res) <- paste(n, '_sum', sep = '')
            oput <- cbind(oput, res)
        }
    }
    
    # compute resultant
    if (summary_type == 'resultant') {
        oput$az <- NA
        oput$mag <- NA
        
        for (i in 1:nrow(oput)) {
            x_sub <- x[month_num == oput$month[i], ]
            res <- calc_resultant (flux = x_sub[, flux_col],
                                          record_duration = x_sub[, record_duration_col],
                                          wind_azimuth = x_sub[, wind_azimuth_col] )
            oput$az[i] <- res$az
            oput$mag[i] <- res$mag
        }
    }
    
    return (oput)
}  



 