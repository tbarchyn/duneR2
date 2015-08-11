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

# PURPOSE: compute some classical style indices of dune mobility.

ecmwf_mobindex <- function (flux, soil_moisture, modifications) {
    # Function to compute a mobility index relating sediment flux to surface
    # soil moisture. Nearly all classical aridity indices are proxies of
    # soil moisture (which in itself it a proxy for vegetation growth). By
    # directly using modeled soil moisture, we avoid much of the issues
    # with aridity metrics (at the expense of comparability). This function
    # only really makes sense on the yearly timescale.
    
    # flux = resultant sediment flux at yearly timescale
    # soil_moisture = a metric of soil moisture
    # modifications = specific modifications (see code)
    
    if (missing (modifications)) {
        modifications <- 'vanilla'
    }
    
    # vanilla mobindex
    if (modifications == 'vanilla') {
        # ratio: higher values indicate more mobility, lower values indicate
        # more vegetation growth.
        mob_index <- flux / soil_moisture
    }

    # vanilla anomaly mobindex: just reports anomalies from the series mean
    if (modifications == 'vanilla_anomaly') {
        # ratio: higher values indicate more mobility, lower values indicate
        # more vegetation growth.
        mob_index <- flux / soil_moisture
        mob_index <- mob_index - mean(mob_index)        # this will fail with NAs
                                                        # which likely isn't a bad thing
    }
    return (mob_index)
}









