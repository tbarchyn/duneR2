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

# PURPOSE: this contains wrapper code for the filter.exe program to
# perform large diameter filter calculations on oversized rasters.

filterR <- function (input, radius, code, output, nontoxic_frac) {
    # Function to run the filter.exe program from R using external GIS
    # files.
    
    # input = input filename
    # radius = the radius of the filter circle in cells
    # code = the operation code to perform:
    #        m = mean
    #        s = sum
    #        f = floor (minimum)
    #        c = ceiling (maximum)
    # output = the output file
    # nontoxic_frac = the minimum acceptable nontoxic fraction (0-1)
    
    # set location to exe file (this will need to be fixed in the future!)
    filter_exe <- 'C://Users//tom//Dropbox//R//duneR//filter//filter.exe'
    
    # set temporary file names (these will be stored in the present working directory)
    asc_input <- paste (input, '.asc', sep = '')
    asc_output <- paste (output, '.asc', sep = '')

    # translate the format of the input file to an arc ascii file
    r <- raster (input)
    r_proj <- projection (r)                
    writeRaster (r, asc_input, overwrite = TRUE)
    
    # assemble callstring for calling filter_exe
    
    callstring <- paste (filter_exe, asc_input, radius, code, asc_output, nontoxic_frac)
    shell (callstring)
    
    r <- raster (asc_output, crs = r_proj)
    writeRaster (r, output, overwrite = TRUE)
    unlink (asc_input)
    unlink (asc_output)
}

calc_residual <- function (raw_topo, radius, nontoxic_frac) {
    # Function to calculate residuals from smoothed minimum surfaces. This automatically
    # creates file names and deals with files on disk. In the future, this will deal
    # with raster objects, but for now, just tiff files on disk.
    
    # raw_topo = raw topography
    # radius = the radius in cells
    # nontoxic_frac = the minimum acceptable nontoxic fraction (0-1)

    # read and write the raw_topo file so the origins match perfectly (there is a rounding
    # issue here with the files that occurs with gdal origins).
    rawt <- raster (raw_topo)
    rawt_proj <- projection (rawt)
    temp_raw_write <- paste (raw_topo, '.asc', sep = '')
    writeRaster (rawt, temp_raw_write, overwrite = TRUE)
    rawt <- raster (temp_raw_write, crs = rawt_proj)
    writeRaster (rawt, raw_topo, overwrite = TRUE)
    unlink (temp_raw_write)
    
    # set up filenames
    min_filename <- paste(radius, '_min.tif', sep = '')
    smin_filename <- paste(radius, '_smin.tif', sep = '')
    rsmin_filename <- paste(radius, '_rsmin.tif', sep = '')

    # compute minimum filter
    filterR (raw_topo, radius, 'f', min_filename, nontoxic_frac)

    # compute smoothed mean filter
    filterR (min_filename, radius, 'm', smin_filename, nontoxic_frac)
    
    # compute residual
    topo <- raster (raw_topo)
    basement <- raster (smin_filename)
    residual <- topo - basement
    writeRaster (residual, rsmin_filename, overwrite = TRUE)
}











