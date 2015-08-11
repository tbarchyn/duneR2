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


# PURPOSE: this software is a collection of relevant R functions that can be strung
# together to do interesting things with data related to dunes. Possibly make it an
# proper R package at some point far far in the future.

# import head: just import this file, and it imports the other files

library (zoo)
library (xts)
library (raster)


source ('./degrib.R')
source ('./utilities.R')
source ('./flux.R')
source ('./filter.R')
source ('./summarize.R')
source ('./mob_indices.R')
