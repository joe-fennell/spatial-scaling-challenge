###########################################################################
###                                                                     ###
####                     SPATIAL SCALING CHALLENGE                     ####
#####                                                                 #####
###########################################################################
# # ------------------------------------------------------------------- # #
# #                             ORGANIZED BY                            # #
# #                                                                     # #
# #  COST ACTION CA17134 "Optical synergies for spatiotemporal SENsing  # #
# #            of Scalable ECOphysiological traits" (SENSECO)           # #
# #                        https://www.senseco.eu                       # #
# #                                                                     # #
# #                            Working Group 1                          # #
# # Closing the scaling gap: from leaf measurements to satellite images # #
# #        https://www.senseco.eu/working-groups/wg1-scaling-gap        # #
# #                                                                     # #
# #     * Dr Javier Pacheco-Labrador (vice-leader)                      # #
# #             Max Planck Institute for Biogeochemistry, Germany       # #
# #     * Dr Ma. Pilar Cendrero-Mateo (leader)                          # #
# #             University of Valencia, Spain                           # #
# #     * Dr Shari Van Wittenberghe, (vice-leader)                      # #
# #             University of Valencia, Spain                           # #
# #     * Dr Gerbrand Koren                                             # #
# #             Utrecht University, Netherlands                         # #
# #     * Prof. Zbynek Malenovský                                       # #
# #             University of Bonn, Germany                             # #
# #                                                                     # #
# # ------------------------------------------------------------------- # #
# #                           CONTACT DETAILS                           # #
# #                                                                     # #
# #   Please, contact us via <scalingchallenge@gmail.com> for any       # #
# #   question or trouble found with the code or the data provided for  # #
# #   the Spatial Scaling Challenge.                                    # #
# #                                                                     # #
# # ------------------------------------------------------------------- # #
# #                      DESCRIPTION AND DISCLAIMER                     # #
# #                                                                     # #
# #   This script opens the netCDF4 and CSV files provided together     # #
# #   this code containing the datasets necessary to participate in     # #
# #   the Spatial Scaling Challenge organized by the COST Action        # #
# #   SENSECO. It also allows exporting the participant's results to    # #
# #   the only standardized format that will be accepted for submission # #
# #                                                                     # #
# #   This program is free software: you can redistribute it and/or     # #
# #   modify it under the terms of the GNU General Public License as    # #
# #   published by the Free Software Foundation, either version 3 of    # #
# #   the License, or any later version.                                # #
# #                                                                     # #
# #   This program is distributed in the hope that it will be useful,   # #
# #   but WITHOUT ANY WARRANTY; without even the implied warranty of    # #
# #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     # #
# #   GNU General Public License for more details.                      # #
# #                                                                     # #
# #   You should have received a copy of the GNU General Public License # #
# #   along with this program. If not, see                              # #
# #   <http://www.gnu.org/licenses/>.                                   # #
# #                                                                     # #
# #   Meteorological data from Majadas de Tiétar experimental station   # #
# #   were provided by the Max Planck Institute for Biogeochemistry     # #
# #   (Germany) and Fundación CEAM (Spain)                              # #
# #                                                                     # #
# # ------------------------------------------------------------------- # #
# #                             CODE AUTHORS                            # #
# #                                                                     # #
# #     * Dr Enrico Tomelleri                                           # #
# #             Free University of Bolzano, Germany                     # #
# #     * Dr Javier Pacheco-Labrador                                    # #
# #             Max Planck Institute for Biogeochemistry, Germany       # #
# #                                                                     # #
# # ------------------------------------------------------------------- # #
# # ------------------------------------------------------------------- # #

## 1) Initialize ------------------------------------------------------ # #
rm(list = ls())

# # load libraries
if(!require(ncdf4)) install.packages("ncdf4")
if(!require(raster)) install.packages("raster")
if(!require(rgdal)) install.packages("rgdal")

# # Options
do_plots <- TRUE

# # Define paths
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

ori_SSCdata <- './1_SSC_data/'
SSC_R <- paste0(ori_SSCdata, 'Airborne_HDRF.nc')
SSC_F <- paste0(ori_SSCdata, 'Airborne_F.nc')
SSC_LST <- paste0(ori_SSCdata, 'Airborne_LST.nc')
SSC_field_plots <- paste0(ori_SSCdata, 'FieldData.csv')
SSC_field_time_series <- paste0(ori_SSCdata, 'FieldData_hh.csv')
ori_SSCresults <- './3_SSC_results/'
out_netcdf <- paste0(ori_SSCresults, 'maps_estimates.nc')

## 2) Spatial Scaling Challenge. Data import -------------------------- # #
# # Coordinates, common to all imagery.

# open a netCDF file
ncin <- nc_open(SSC_R)
spat_ref <- ncatt_get(ncin, 0, "spat_ref")$value
pixel_size <- ncatt_get(ncin, 0, "pixel_size")$value#in meters

coords_x <- ncvar_get(ncin,"coords_x") + pixel_size/2#from top-right corner to center of the pixel
coords_y <- ncvar_get(ncin,"coords_y") - pixel_size/2#from top-right corner to center of the pixel
coords_z <- ncvar_get(ncin,"coords_z")
n_rows <- dim(coords_x)[1]
n_cols <- dim(coords_x)[2]

# # Airborne hyperspectral reflectance imagery
R_ <- ncvar_get(ncin, 'var_')#(Reflectance imagery [n_rows, n_cols, n_bands])
R_wvl <- ncvar_get(ncin, 'wvl_center')#(Center wavelength)
R_fwhm <- ncvar_get(ncin, 'wvl_fwhm')#(Full Width Half Maximum)
R_sza = ncvar_get(ncin, 'VZA')#View Zenith Angle
R_svaa = ncvar_get(ncin, 'SVAA')#Sun-View Azimuth Angle
R_description <- ncatt_get(ncin, 0, 'var_description')$value
R_units <- ncatt_get(ncin, 0, 'var_uds')$value
R_wvl_units <- ncatt_get(ncin, 0, 'wvl_uds')$value
n_bands <- length(R_wvl)
nc_close(ncin)

# # Airborne sun-induced chlorophyll fluorescence radiance imagery
ncin <- nc_open(SSC_F)
F_ <- ncvar_get(ncin, 'var_')#(Fluorescence radiance imagery [n_rows, n_cols, n_bands])
F_wvl <- ncvar_get(ncin, 'wvl_center')#(Center wavelength)
F_fwhm <- ncvar_get(ncin, 'wvl_fwhm')#(Full Width Half Maximum)
F_sza = ncvar_get(ncin, 'VZA')#View Zenith Angle
F_svaa = ncvar_get(ncin, 'SVAA')#Sun-View Azimuth Angle
F_description <- ncatt_get(ncin, 0, 'var_description')
F_units <- ncatt_get(ncin, 0, 'var_uds')
F_wvl_units <- ncatt_get(ncin, 0, 'wvl_uds')
nc_close(ncin)

# # Airborne land surface temperature imagery
ncin <- nc_open(SSC_LST)
LST_ <- ncvar_get(ncin, 'var_')#(Land surface temperature imagery [n_rows, n_cols])
LST_description <- ncatt_get(ncin, 0, 'var_description')
LST_sza = ncvar_get(ncin, 'VZA')#View Zenith Angle
LST_svaa = ncvar_get(ncin, 'SVAA')#Sun-View Azimuth Angle
LST_units <- ncatt_get(ncin, 0, 'var_uds')
nc_close(ncin)

# # Field data. Spatial sampling in 1 x 1 meter plots.
FP_ <- read.csv(SSC_field_plots)

# # Field data. Meteorological data and time series of NPQ.
FPhh_ <- read.csv(SSC_field_time_series)

if (do_plots == TRUE){
  # # Airborne hyperspectral reflectance imagery
  bsel_ <- which(abs(R_wvl-680) == min(abs(R_wvl-680)))
  ras_plot <- rasterFromXYZ(data.frame(x = as.vector(coords_x),
                        y = as.vector(coords_y),
                        z = as.vector(R_[, , bsel_])))
  pal <- colorRampPalette(c("green2","yellow2"))
  bbox_ <- bbox(ras_plot)
  points_plots <- SpatialPoints(data.frame(x=FP_$Xutm, y=FP_$Yutm), bbox=bbox_)
  points_plots_lab <- c(sprintf("%s", FP_$PlotNum[]))
  points_npq_lab = c("moni-PAM                ", "                 moni-PAM")
  points_npq <- SpatialPoints(
    data.frame(x=c(FPhh_$XNPQwheat, FPhh_$XNPQmaize),
               y=c(FPhh_$YNPQwheat, FPhh_$YNPQmaize)), bbox=bbox_)
  
  png("1_HDRF_680.png", width = 675, height = 700)
  plot(ras_plot, col=pal(7), xlab = "x (m)", ylab = "y (m)")
  plot(points_plots, add = TRUE, pch = 16, col='red')
  text(points_plots, labels = points_plots_lab, pos = 4, offset = 0.7, col='red')
  plot(points_npq, add = TRUE, pch = 16, col='blue')
  text(points_npq, labels = points_npq_lab, pos = 1, offset = .7, col='blue')
  text(SpatialPoints(data.frame(
    x=coords_x[1,1]+10, y=coords_y[1,1]-1), bbox=bbox_),
    labels = "Field 1 (Triticum aestivum)", pos = 4, offset = .7)
  text(SpatialPoints(data.frame(
    x=coords_x[1,1]+65, y=coords_y[1,1]-1), bbox=bbox_),
    labels = "Field 2 (Zea mays)", pos = 4, offset = .7)
   title(sprintf('HDRF @ 680 nm (%s)', R_units))
   dev.off()

   # # Airborne sun-induced chlorophyll fluorescence radiance imagery
   ras_plot <- rasterFromXYZ(data.frame(x = as.vector(coords_x),
                                        y = as.vector(coords_y),
                                        z = as.vector(F_[, , 2])))
   png("2_F_760.png", width = 675, height = 700)
   plot(ras_plot, col=pal(7), xlab = "x (m)", ylab = "y (m)")
   plot(points_plots, add = TRUE, pch = 16, col='red')
   text(points_plots, labels = points_plots_lab, pos = 4, offset = 0.7, col='red')
   plot(points_npq, add = TRUE, pch = 16, col='blue')
   text(points_npq, labels = points_npq_lab, pos = 1, offset = .7, col='blue')
   text(SpatialPoints(data.frame(
     x=coords_x[1,1]+10, y=coords_y[1,1]-1), bbox=bbox_),
     labels = "Field 1 (Triticum aestivum)", pos = 4, offset = .7)
   text(SpatialPoints(data.frame(
     x=coords_x[1,1]+65, y=coords_y[1,1]-1), bbox=bbox_),
     labels = "Field 2 (Zea mays)", pos = 4, offset = .7)
   title(sprintf('F @ 760 nm (%s)', F_units))
   dev.off()

   # # Airborne Land Surface Temperature
   ras_plot <- rasterFromXYZ(data.frame(x = as.vector(coords_x),
                                        y = as.vector(coords_y),
                                        z = as.vector(LST_)))
   png("3_LST.png", width = 675, height = 700)
   plot(ras_plot, col=pal(7), xlab = "x (m)", ylab = "y (m)")
   plot(points_plots, add = TRUE, pch = 16, col='red')
   text(points_plots, labels = points_plots_lab, pos = 4, offset = 0.7, col='red')
   plot(points_npq, add = TRUE, pch = 16, col='blue')
   text(points_npq, labels = points_npq_lab, pos = 1, offset = .7, col='blue')
   text(SpatialPoints(data.frame(
     x=coords_x[1,1]+10, y=coords_y[1,1]-1), bbox=bbox_),
     labels = "Field 1 (Triticum aestivum)", pos = 4, offset = .7)
   text(SpatialPoints(data.frame(
     x=coords_x[1,1]+65, y=coords_y[1,1]-1), bbox=bbox_),
     labels = "Field 2 (Zea mays)", pos = 4, offset = .7)
   title(sprintf('LST (%s)', LST_units))
   dev.off()
}

## 3) Brief example of how to link field and imagery data ------------- # #
# Generate indices for each crop
Iw <- (FP_$Crop == 'Wheat')
Im <- !Iw

# # Access to the first layer of the [row, col, band] imagery matrix
if (do_plots == TRUE){png("4_F687_vs_LAI.png", width = 675, height = 700)}
ind_ <- cbind(FP_$Ypix, FP_$Xpix,matrix(1, nrow = length(FP_$Xpix), ncol = 1))
plot(F_[ind_[Iw, ]], FP_$LAI[Iw], col='red', ylim=c(2,6), xlim=c(.9,1.4),
     xlab = "F @ 687 nm",  ylab="LAI (m2 m-2)")
points(F_[ind_[Im, ]], FP_$LAI[Im], col='blue')
if (do_plots == TRUE){dev.off()}


# # Access to the second layer of the [row, col, band] imagery matrix
if (do_plots == TRUE){png("5_F760_vs_LAI.png", width = 675, height = 700)}
ind_ <- cbind(FP_$Ypix, FP_$Xpix,matrix(2, nrow = length(FP_$Xpix), ncol = 1))
plot(F_[ind_[Iw, ]], FP_$LAI[Iw], col='red', ylim=c(2,6), xlim=c(1.1,2),
     xlab= "F @ 760 nm", ylab="LAI (m2 m-2)")
points(F_[ind_[Im, ]], FP_$LAI[Im], col='blue')
if (do_plots == TRUE){dev.off()}

# # Access to the n-layer layer of the [row, col, band] imagery matrix
if (do_plots == TRUE){png("6_HDRF680_vs_LAI.png", width = 675, height = 700)}
bsel_ <- which(abs(R_wvl-680) == min(abs(R_wvl-680)))
ind_ <- cbind(FP_$Ypix, FP_$Xpix,matrix(bsel_, nrow = length(FP_$Xpix), ncol = 1))
plot(R_[ind_[Iw, ]], FP_$LAI[Iw], col='red', ylim=c(2,6), xlim=c(0.02,0.06),
     xlab= "HDRF @ 680 nm", ylab="LAI (m2 m-2)")
points(R_[ind_[Im, ]], FP_$LAI[Im], col='blue')
if (do_plots == TRUE){dev.off()}

###########################################################################
# # YOUR WORK STARTS HERE!

## 4) Your turn. Solve the challenge! --------------------------------- # #



## 5) Your turn. Describe your methods -------------------------------- # #
# Before submitting your results, you need to document the methods you
# used in a standardized way that will help to summarize the contribution
# of all the participants. Use any of the MS Word (.doc or .docx) or
# OpenOffice Writer (.odt) templates: /2_SSC_templates/SSC_report.xxx
# First, provide your diagnosis of the actual vegetation status. Then, 
# describe the methods you used to estimate each of the variables (and
# uncertainties if you did). Try to be concise and clear. The descriptions
# could be included in the supplementary material of the joint manuscript,
# therefore, take care of English grammar and style.
# Include references if necessary (Author et al., year) and the full
# reference at the end of each section.

# # Here, select the filename matching the extension of your methods' file
out_methods = paste0(ori_SSCresults,'SSC_report.docx')
# out_methods = [ori_SSCresults,'SSC_report.doc'];
# out_methods = [ori_SSCresults,'SSC_report.odt'];

## 6) Your turn. Prepare your results --------------------------------- # #
# # 6.1) Fill up this structure with your personal data
pdata_ <- list()
pdata_$name <- ''
pdata_$middle_name <- ''
pdata_$surname <- ''
pdata_$orcid <- ''
pdata_$institution <- ''
pdata_$department <- ''
pdata_$address <- ''
pdata_$email <- ''
pdata_$positionexperience <- ''
# This is the name that will be used for the compressed zip folder that
# will be generate. Normally, it should just be the surname followed by
# the name. However, if your surname or your name included characters
# that could not be used for a folder name, write a valid version instead.
pdata_$surname_name4filename <- paste0(pdata_$surname,pdata_$name)

# # 6.2) Fill up this structure with your estimates and, if you did it,
# with the estimated uncertainties of your maps. Important, these
# variables must have the same dimensions [n_rows, n_cols] than the imagery
# provided, and in the units decribed below
results <- list()

# # Estimates
results$LAI_est <- matrix(-999., n_rows, n_cols)# Estimated map of leaf area index [m^2 m^-2]
results$Cab_est <- matrix(-999., n_rows, n_cols)# Estimated map of leaf chlorophyll content [ug cm^-2]
results$Vcmax25_est <- matrix(-999., n_rows, n_cols)# Estimated map of maximum carboxylation rate at 25 ?C [umol cm^-2 s^-1]
results$NPQ_est <- matrix(-999., n_rows, n_cols)# Estimated map of maximum carboxylation rate at 25 ?C [umol cm^-2 s^-1]

# # Uncertainties (leave -999. if you did not estimate them)
results$LAI_unc <- matrix(-999., n_rows, n_cols)# Estimated uncertainties of leaf area index [m^2 m^-2]
results$Cab_unc <- matrix(-999., n_rows, n_cols)# Estimated uncertainties of leaf chlorophyll content [ug cm^-2]
results$Vcmax25_unc <- matrix(-999., n_rows, n_cols)# Estimated uncertainties of maximum carboxylation rate at 25 ?C [umol cm^-2 s^-1]
results$NPQ_unc <- matrix(-999., n_rows, n_cols)# Estimated uncertainties of maximum carboxylation rate at 25 ?C [umol cm^-2 s^-1]

# # Stress map range between 0 for minimum stress and 1 for maximum stress
results$stress_est <- matrix(-999., n_rows, n_cols)# % Estimated map of stress [-]

# # PERFECT!!
# Now wait for the zip file to be generated and send it by email to the
# email provided in the contact section <scalingchallenge@gmail.com>


###########################################################################
# # YOUR WORK ALMOST FINISHES HERE!

## 7) Spatial Scaling Challenge. Prepare standard output files -------- # #
# # Check if the methods section is included the results folder
if (file.exists(out_methods) == FALSE || pdata_$surname_name4filename == ''){
    warning('The results are not still ready. The zip file will not be produced')
} else {
# # Remove first the output file if it exists
    if (file.exists(out_netcdf) == TRUE){
      file.remove(out_netcdf)
    }

# # Add results
   xdim <- ncdim_def("x","m_east", coords_x[1, ]) 
   ydim <- ncdim_def("y","m_north",coords_y[1, ])
   LAI_est <- ncvar_def(name='LAI_est', units='m^2 m^-2', dim=list(xdim, ydim), missval = -99., longname='Estimated map of leaf area index', prec='double')
   Cab_est <- ncvar_def(name='Cab_est', units='ug cm^-2', dim=list(xdim, ydim), missval = -99., longname='Estimated map of leaf chlorophyll', prec='double')
   Vcmax25_est <- ncvar_def(name='Vcmax25_est', units='umol cm^-2 s^-1', dim=list(xdim, ydim), missval = -99.,  longname='Estimated map of maximum carboxylation rate at 25 ?C', prec='double')
   NPQ_est <- ncvar_def(name='NPQ_est', units='umol cm^-2 s^-1', dim=list(xdim, ydim), missval = -99., longname='Estimated map of maximum carboxylation rate at 25 ?C', prec='double')
   LAI_unc <- ncvar_def(name='LAI_unc', units='m^2 m^-2', dim=list(xdim, ydim),  missval = -99., longname='Estimated uncertainties of leaf area index', prec='double')
   Cab_unc <- ncvar_def(name='Cab_unc', units='ug cm^-2', dim=list(xdim, ydim), missval = -99., longname='Estimated uncertainties of leaf chlorophyll content', prec='double')
   Vcmax25_unc <- ncvar_def(name='Vcmax25_unc', units='umol cm^-2 s^-1',  dim=list(xdim, ydim), missval = -99., longname='Estimated uncertainties of maximum carboxylation rate at 25 ?C', prec='double')
   NPQ_unc <- ncvar_def(name='NPQ_unc', units='umol cm^-2 s^-1', dim=list(xdim, ydim),  missval = -99., longname='Estimated uncertainties of maximum carboxylation rate at 25 ?C', prec='double')
   stress_est <- ncvar_def(name='stress_est', units='-', dim=list(xdim, ydim), missval = -99., longname='Estimated map of stress',prec='double')
   Vars  <- list(LAI_est, Cab_est, Vcmax25_est, NPQ_est, LAI_unc,
                 Cab_unc, Vcmax25_unc, NPQ_unc, stress_est)

# # Create a new netcdf file and insert variables.
    con <- nc_create(out_netcdf, Vars, force_v4 = TRUE, verbose = TRUE)
    ncvar_put(con, LAI_est, results$LAI_est, start=c(1, 1))
    ncvar_put(con, Cab_est, results$Cab_est)
    ncvar_put(con, Vcmax25_est, results$NPQ_est)
    ncvar_put(con, NPQ_est, results$NPQ_est)
    ncvar_put(con, LAI_unc, results$LAI_unc)
    ncvar_put(con, Cab_unc, results$Cab_unc)
    ncvar_put(con, Vcmax25_unc, results$Vcmax25_unc)
    ncvar_put(con, NPQ_unc, results$NPQ_unc)
    ncvar_put(con, stress_est, results$stress_est)
    
# # Add personal data
    ncatt_put(con,0,"name",pdata_$name)
    ncatt_put(con,0,"middle_name",pdata_$middle_name)
    ncatt_put(con,0,"surname",pdata_$surname)
    ncatt_put(con,0,"orcid",pdata_$orcid)
    ncatt_put(con,0,"institution",pdata_$institution)
    ncatt_put(con,0,"department",pdata_$department)
    ncatt_put(con,0,"address",pdata_$address)
    ncatt_put(con,0,"email",pdata_$email)
    ncatt_put(con,0,"surname_name4filename",pdata_$surname_name4filename)
    ncatt_put(con,0,"positionexperience",pdata_$positionexperience)
    ncatt_put(con,0,"Scrip_SSC",'R')
    
# # close the file for having it written to disk
    nc_close(con)    
    
# # Compress the files to be sent to SENSECO Working Group 1
    zip_fname <- paste0(pdata_$surname_name4filename, '_4sub.zip')
    zip(zipfile=paste0(ori_SSCresults, zip_fname), files=c(out_netcdf, out_methods))

# # Final instructions
    print('Congratulations!! Your results have been stored in the following zip file:\n')
    print(zip_fname)
    print('Send this file via email to the Spatial Scaling Challenge email adress')
    print('<scalingchallenge@gmail.com>')
    print('Thanks a lot for your participation!')
    print('Looking forward to learn from the community and share the manuscript draft any soon!')
    print('COST Action SENSECO Working Group')
}
