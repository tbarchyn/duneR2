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

# PURPOSE: these functions are generic worker functions.

uvw_to_scalar <- function (u, v, w, convention) {
    # Function to convert u, v, w wind components to wind direction, speed,
    # and altitude. The wind direction is the direction that the wind is COMING
    # FROM, and the speed is the 3 dimensional speed. If the w component
    # is missing, of course the wind speed will be only 2d. Set the w
    # component to 0 if using u, v, w data and you desire only 2d wind speed.
    
    # Returns a dataframe with wspd, wdir, and alt columns for wind speed, wind
    # direction, and wind altitude.
    
    # There are several conventions for u, v, w data.  See notes
    # below for what the u, v, w wind components mean when using each respective
    # convention.
    
    # convention: '81000_sonic', 'climate'
    
    if (missing (w)) {
        w <- 0.0            # allow missing w, by setting to 0.0
    }
    
    # force the inputs to be vectors, this allows random timeseries classes
    # to be given as inputs and results in a regularly formatted data.frame output.
    u <- as.vector (u)
    v <- as.vector (v)
    w <- as.vector (w)
    
    if (convention == '81000_sonic') {
        # The conventions for the 81000 sonic anemometer are:
        # u = the east-west wind speed, wind from the east is positive u,
        #     wind from the west is negative u.
        # v = the north-south wind speed, wind from the north is positive v,
        #     wind from the south is negative v.
        # w = the up-down wind speed, wind from below is positive, wind from
        #     above is negative.        # calculate wind direction, and deal with cases where azimuth is out of 0-360.0
        
        # calculate wind direction, and deal with cases where azimuth is out of 0-360.0
        wind_dir <- 90.0 - (atan2(v, u) * (180.0 / pi))
        wind_dir[!is.na(wind_dir) & wind_dir >= 360.0] <- wind_dir[!is.na(wind_dir) & wind_dir >= 360.0] - 360.0
        wind_dir[!is.na(wind_dir) & wind_dir < 0.0] <- wind_dir[!is.na(wind_dir) & wind_dir < 0.0] + 360.0
    
        if (length(w[w != 0.0]) > 0) {
            # calculate 3D wind speed
            wind_speed <- sqrt (u^2 + v^2 + w^2)
            
            # calculate altitude, positive when wind is coming from below
            altitude <- atan ((w / sqrt(u^2 + v^2))) * (180.0 / pi)
        } else {
            # calculate 2D wind speed
            wind_speed <- sqrt (u^2 + v^2)
            
            # set altitude
            altitude <- rep (0.0, length (wind_speed))
        }
    } else if (convention == 'climate') {
        # The conventions for most climate data are:
        # u = the east-west wind speed, wind from the west is positive u,
        #     wind from the east is negative u.
        # v = the north-south wind speed, wind from the south is positive v,
        #     wind from the north is negative v.
        # w = the up-down wind speed, wind from below is positive, wind from
        #     above is negative.               
        
        # calculate wind direction, and deal with cases where azimuth is out of 0-360.0
        wind_dir <- atan2((-1.0 * u), (-1.0 * v)) * (180.0 / pi)
        wind_dir[!is.na(wind_dir) & wind_dir >= 360.0] <- wind_dir[!is.na(wind_dir) & wind_dir >= 360.0] - 360.0
        wind_dir[!is.na(wind_dir) & wind_dir < 0.0] <- wind_dir[!is.na(wind_dir) & wind_dir < 0.0] + 360.0
        
        if (length(w[w != 0.0]) > 0) {
            # calculate 3D wind speed
            wind_speed <- sqrt (u^2 + v^2 + w^2)
            
            # calculate altitude, positive when wind is coming from below
            altitude <- atan ((w / sqrt(u^2 + v^2))) * (180.0 / pi)
        } else {
            # calculate 2D wind speed
            wind_speed <- sqrt (u^2 + v^2)
            
            # set altitude
            altitude <- rep (0.0, length (wind_speed))
        }
    } else {
        print ('ERROR: undefined convention for wind components')
    }
    
    output <- data.frame (wdir = wind_dir, wspd = wind_speed, alt = altitude)
    return (output)
}

calc_u_star_simple <- function (wspd, anemo_elev, z0) {
    # Function to calculate u_star with one elevation of anemometer and a z0
    # value. Note that u_star is best calculated with a wind profile, rather
    # than just one measurement.
    
    # wspd = the wind speed measurement (m/s)
    # anemo_elev = anemometer elevation above the surface (m)
    # z0 = the aerodynamic roughness (z0) value assumed (m)
    
    wspd <- as.vector (wspd)
    anemo_elev <- as.vector (anemo_elev)
    z0 <- as.vector (z0)
    
    von_karman <- 0.41
    u_star <- wspd / ((1.0 / von_karman) * log(anemo_elev / z0))
    
    return (u_star)
}

calc_wspd <- function (u_star, elev, z0) {
    # Function to calculate a wind speed at a specific elevation from a given
    # u_star with the 'Law of the Wall'
    
    # u_star = the friction velocity (m/s)
    # elev = the elevation to return new wind speed (m)
    # z0 = the aerodynamic roughness (z0) value assumed (m)
    
    u_star <- as.vector (u_star)
    elev <- as.vector (elev)
    z0 <- as.vector (z0)

    von_karman <- 0.41
    wspd <- (u_star / von_karman) * log(elev / z0)
    
    return (wspd)
}

calc_mean_wind_dir <- function (u, v, convention) {
    # Function to calculate mean wind direction
    # u = u wind speeds
    # v = v wind speeds
    mean_u <- mean(u, na.rm = T)
    mean_v <- mean(v, na.rm = T)
    direction <- uvw_to_scalar (u = mean_u, v = mean_v, w = 0.0, convention)$wdir
    return (direction)
}

calc_air_density <- function (pressure, kelvin_temp, xts_interpolate) {
    # Function to calculate the air density from pressure and air temperature.
    
    # pressure = the pressure in Pa
    # kelvin_temp = the temperature in K
    # xts_interpolate = interpolate the NAs in the timeseries (as is necessary with
    #                   ecmwf data because temps are only every 6 hours, where the
    #                   pressures are every 3 hours
    
    pressure <- as.vector (pressure)
    
    if (missing (xts_interpolate)) {
        xts_interpolate <- FALSE
    }
    
    if (length(kelvin_temp[kelvin_temp < 200]) > 0) {
        print ('WARNING: calc_air_density possibly received temperature in C, not K')
    }
    
    # xts interpolate
    if (xts_interpolate) {
        # call the xts_interpolate wrapper function
        kelvin_temp <- xts_interpolate (kelvin_temp)
    } else {
        # else, just ensure it is a vector, not an xts object
        kelvin_temp <- as.vector (kelvin_temp)
    }
    
    R_air <- 287.058
    
    rho_air <- pressure / (R_air * kelvin_temp)
    return (rho_air)
}

xts_interpolate <- function (xts_in) {
    # Wrapper function for the na.approx function in zoo. The problem is that if the last
    # value in a timeseries is a NA, na.approx returns a xts object that is 1 record shorter
    # and wreaks havoc downstream. This function should be used sparingly and carefully.
    
    original_len <- length(xts_in)
    xts_in <- na.approx(xts_in)
    vec <- as.vector(xts_in)
    if (original_len - length(xts_in) == 1) {
        # the na.approx function can fail to interpolate the last point, and produce
        # a timeseries which lacks the same number of points as the input values
        # to rectify this, we append 1 value to the end.
        vec <- c(vec, vec[length(vec)])
    }
    return (vec)
}

xts_index <- function (xts_in) {
    # Convenience function to create character strings of the index of xts
    # objects that have proper millisecond values. Often the time index will
    # round improperly during conversion to an xts object, and thus create
    # issues for merging, and it will print incorrectly with strptime due to
    # decimal place truncation (not rounding!). This just outputs a pseudo standard
    # POSIX character string that will read nicely into excel or backconvert
    # straightforwardly with as.POSIX functions.
    
    # xts_in = the input xts object
    
    char_index <- strftime (index (xts_in), '%Y-%m-%d %H:%M:')
    sec_index <- strftime (index (xts_in), '%OS6')
    sec_index <- as.numeric (sec_index)
    sec_index <- round(sec_index, 1)
    char_index <- paste(char_index, as.character(sec_index), sep = '')
    return (char_index)
}

reload_duneR <- function () {
    # Convenience function to reload duneR during development (remove later!)
    original_dir <- getwd()
    setwd ('C:/Users/tom/Dropbox/R/duneR2')
    source ('duneR.R')
    setwd (original_dir)
}


