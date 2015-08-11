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

# PURPOSE: these functions are for calculating sediment flux and doing
# things with those calculated sediment fluxes.

calc_u_star_t <- function (..., method) {
    # Wrapper function for calculating the threshold shear velocity for
    # sand transport. The method is necessary to call a specific method
    # for calculation. Don't see the method you want to use? Please
    # go ahead and add it.
    
    # ... = the parameters passed to the method
    # method = the method desired

    if (method == 'shao_lu') {
        # use the shao_lu method for calculating the threshold shear velocity
        # in m/s. Refer to calc_u_star_shao_lu function for details.
        u_star_t <- calc_u_star_t_shao_lu (...)
    
    } else {
        print ('ERROR: undefined threshold calculation method')
    }
    return (u_star_t)
}

calc_flux <- function (..., method) {
    # Wrapper function for calculating the aeolian sediment flux. The
    # subset methods are called here, but defined elsewhere.
    # Don't see the method you want to use? Please go ahead and add it.

    # ... = the parameters passed to the method
    # method = the method desired

    if (method == 'corrected_white') {
        # use the corrected_white method to calculate flux in kg/s per
        # cross-wind meter, refer to calc_flux_corrected_white.
        flux <- calc_flux_corrected_white (...)
    
    } else if (method == 'corrected_white_bulk') {
        # use the corrected_white method to calculate flux in kg/s per
        # cross-wind meter, refer to calc_flux_corrected_white for 
        # details. Then, take the output and calculate the bulk flux
        # in m^3/s per cross-wind meter. This is a convenience function
        # that assumes bulk density of 1500 kg/m^3 (note this is different
        # from mineral density of most aeolian sediment, the difference is
        # the porosity of the sediment).
        
        bulk_density <- 1500.0          # bulk density of in place aeolian
                                        # sediment (kg / m^3)
        flux_mass <- calc_flux_corrected_white (...)
        flux <- flux_mass / bulk_density

    } else {
        print ('ERROR: undefined flux calculation method')
    }
    return (flux)
}

calc_u_star_t_shao_lu <- function (d, rho_sed, rho_air) {
    # Function to calculate the u_star_t with the Shao and Lu (2000) method:

    # Shao, Y., Lu, H., 2000, A simple expression for wind erosion threshold
    # friction velocity, Journal of Geophysical Research 105,
    # 22,437-22,443.
    
    # d = the median sediment diameter (m)
    # rho_sed = the density of the particles (kg / m^3)
    # rho_air = the density of the air (kg / m^3)
    
    # set empirical constants 
    A_n <- 0.0123
    gma <- 0.0003
    g <- 9.81           # acceleration of gravity (m / s^2)
    
    d <- as.vector (d)
    rho_sed <- as.vector (rho_sed)
    rho_air <- as.vector (rho_air)
    
    u_star_t <- sqrt (A_n * ((rho_sed * g * d) + (gma / (rho_air * d))))
    return (u_star_t)
}

calc_flux_corrected_white <- function (u_star, u_star_t, rho_air, force_threshold) {
    # Function to calculate the sediment flux with the corrected White (1979)
    # model. Returns the flux in kg/s for 1 m crosswind. 
    
    # This equation was originally published in: 
        
    # White, B. R., 1979., Soil transport by winds on Mars, Journal of
    # Geophysical Research, 84, 4643–4651.
    
    # The model was subsequently corrected in:
    
    # Namikas, S., Sherman, D.J., 1997, Predicting aeolian sand transport:
    # revisiting the White model, Earth Surface Processes and Landforms, 22,
    # 601-604.
    
    # u_star = friction velocity (m / s)
    # u_star_t = friction velocity fluid threshold for sediment transport (m / s)
    # rho_air = air density (kg / m^3)
    # force_threshold = boolean to enforce the threshold. If true, it makes any
    #                   records where u_star < u_star_t forced to 0.0.
    
    u_star <- as.vector (u_star)
    u_star_t <- as.vector (u_star_t)
    rho_air <- as.vector (rho_air)
    
    if (missing (force_threshold)) {
        force_threshold <- TRUE
    }
    
    g <- 9.81           # accelleration of gravity in m / s^2
    
    flux <- ( 2.61 * (rho_air / g) * u_star^3 * (1.0 - (u_star_t / u_star)) *
            ((1.0 + (u_star_t / u_star))^2) )
    
    # force the threshold if desired
    if (force_threshold) {
        flux[u_star < u_star_t] <- 0.0
    }
    
    # do not return negative fluxes (this is not physically reasonable)
    flux[flux < 0.0] <- 0.0
    
    return (flux)
}    
    
calc_resultant <- function (flux, record_duration, wind_azimuth) {
    # Function to calculate the resultant sediment flux from a series of flux
    # estimates, and the wind direction for each flux estimate. This is simply
    # the vector sum. The resultant is the downwind direction of flux and is
    # composed of a downwind vector and magnitude. Note that the magnitude will
    # be less than a straight sum of sediment flux, and certain dune mobility
    # indices will need to take care what is used (e.g., for a seasonally bimodal
    # wind regime, you likely don't want to use a vector sum flux).
    
    # flux = the average flux rate for each measurement (/s)
    # record_duration = the duration of record (s)
    # wind azimuth = the direction the wind is coming from (degrees)
    
    # calculate the u and v components of the flux with climate conventions
    u_flux <- -1.0 * sin(wind_azimuth * pi / 180.0) * flux * record_duration
    v_flux <- -1.0 * cos(wind_azimuth * pi / 180.0) * flux * record_duration
    
    # sum u and v flux components, this will fail with NAs, which is probably good
    # so the user has to figure out what to do with NAs.
    u_flux_sum <- sum (u_flux)
    v_flux_sum <- sum (v_flux)
    
    # calculate resultant flux angle (this is the downwind angle)
    direction <- atan2(u_flux_sum, v_flux_sum) * (180.0 / pi)
    if (direction < 0.0) {
        direction <- direction + 360.0
    }
    magnitude <- sqrt (u_flux_sum^2 + v_flux_sum^2)
    output <- data.frame (az = direction, mag = magnitude)
    return (output)
}

    

    