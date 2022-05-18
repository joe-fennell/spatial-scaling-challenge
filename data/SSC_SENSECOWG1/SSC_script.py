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
# #     * Dr Gerbrand Koren                                             # #
# #             Utrecht University, Netherlands                         # #
# #     * Dr Javier Pacheco-Labrador                                    # #
# #             Max Planck Institute for Biogeochemistry, Germany       # #
# #                                                                     # #
# # ------------------------------------------------------------------- # #
# # ------------------------------------------------------------------- # #

#%% 1) Initialize

# # Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import netCDF4 as nc
import pandas as pd
from os.path import exists
from os import remove
import warnings
from zipfile import ZipFile

# # Options
do_plots = True

# # Define paths
ori_SSCdata = '1_SSC_data//'
SSC_R = ori_SSCdata + 'Airborne_HDRF.nc'
SSC_F = ori_SSCdata + 'Airborne_F.nc'
SSC_LST = ori_SSCdata + 'Airborne_LST.nc'
SSC_field_plots = ori_SSCdata + 'FieldData.csv'
SSC_field_time_series = ori_SSCdata + 'FieldData_hh.csv'
ori_SSCresults = '3_SSC_results//'
out_netcdf = ori_SSCresults + 'maps_estimates.nc'

#%% 2) Spatial Scaling Challenge. Data import
# # Coordinates, common to all imagery. 
h = nc.Dataset(SSC_R)
spat_ref = h.getncattr('spat_ref')
pixel_size = h.getncattr('pixel_size') #in meters
coords_x =  h.variables['coords_x'][:] + pixel_size/2 #from top-right corner to center of the pixel
coords_y =  h.variables['coords_y'][:] - pixel_size/2 #from top-right corner to center of the pixel
coords_z =  h.variables['coords_z'][:]
[n_rows, n_cols] = np.shape(coords_x)
h.close()

# # Airborne hyperspectral reflectance imagery
h = nc.Dataset(SSC_R)
R_ = np.transpose(h.variables['var_'][:],(1, 2, 0)) #(Reflectance imagery [n_rows, n_cols, n_bands])
R_wvl = h.variables['wvl_center'][:] #(Center wavelength)
R_fwhm = h.variables['wvl_fwhm'][:] #(Full Width Half Maximum)
R_sza = h.variables['VZA'] #(View Zenith Angle)
R_svaa = h.variables['SVAA'] #(Sun-View Azimuth Angle)
R_description = h.getncattr('var_description')
R_units = h.getncattr('var_uds')
R_wvl_units = h.getncattr('wvl_uds')
n_bands = len(R_wvl)
h.close()

# # Airborne sun-induced clorophyll fluorescence radiance imagery
g = nc.Dataset(SSC_F)
F_ = np.transpose(g.variables['var_'][:],(1, 2, 0)) #(Fluorescence radiance imagery [n_rows, n_cols, n_bands])
F_wvl = g.variables['wvl_center'][:] #(Center wavelength)
F_fwhm = g.variables['wvl_fwhm'][:] #(Full Width Half Maximum)
F_sza = g.variables['VZA'] #(View Zenith Angle)
F_svaa = g.variables['SVAA'] #(Sun-View Azimuth Angle)
F_description = g.getncattr('var_description')
F_units = g.getncattr('var_uds')
F_wvl_units = g.getncattr('wvl_uds')
g.close()

# # Airborne land surface temperature imagery
q = nc.Dataset(SSC_LST)
LST_ = np.transpose(q.variables['var_'][:],(1, 2, 0))[:, :, 0] #(Land surface temperature imagery, [n_rows, n_cols])
LST_sza = q.variables['VZA'] #(View Zenith Angle)
LST_svaa = q.variables['SVAA'] #(Sun-View Azimuth Angle)
LST_description = q.getncattr('var_description')
LST_units = q.getncattr('var_uds')
q.close()

# # Field data. Spatial sampling in 1 x 1 meter plots.
FP_ = pd.read_csv(SSC_field_plots, delimiter=',')

# # Field data. Meteorological data and time series of NPQ.
FPhh_ = pd.read_csv(SSC_field_time_series, delimiter=',')

if do_plots is True:    
    # # Airborne hyperspectral reflectance imagery
    bsel_ = np.where(abs(R_wvl-680) == min(abs(R_wvl-680)))[0][0]
    fig = plt.figure(figsize =(10.5, 9))
    ax = plt.axes()
    ax.pcolormesh(coords_x, coords_y, R_[:, :, bsel_], shading='auto',
                  cmap=cm.summer)
    plt.plot(FP_.Xutm, FP_.Yutm, 'or', markersize=11 )    
    for i_ in range(len(FP_)):
        plt.text(FP_.Xutm[i_] + 3*pixel_size, FP_.Yutm[i_] + 3*pixel_size,
                 ('%d' % FP_.PlotNum[i_]), color='r', fontsize=15,
                 fontweight='bold')
    plt.plot(FPhh_.XNPQwheat[0], FPhh_.YNPQwheat[0], 'sb', markersize=11)  
    plt.text(FPhh_.XNPQwheat[0] - 20*pixel_size, FPhh_.YNPQwheat[0],
                 'moni-PAM', color='b', fontsize=15, fontweight='bold')
    plt.plot(FPhh_.XNPQmaize[0], FPhh_.YNPQmaize[0], 'sb', markersize=11)  
    plt.text(FPhh_.XNPQmaize[0] + 3*pixel_size, FPhh_.YNPQmaize[0],
                 'moni-PAM', color='b', fontsize=15, fontweight='bold')
    plt.text(coords_x[0, 0] + n_cols/4,coords_y[0, 0] + 5,
        r'Field 1 (${Triticum}$ ${aestivum}$)', horizontalalignment='center',
        fontsize=12)
    plt.text(coords_x[0, 0] + 3*n_cols/4,coords_y[0, 0] + 5,
        r'Field 2 (${Zea}$ ${mays}$)', horizontalalignment='center',
        fontsize=12)
    plt.xlabel(r'$x$ (m)')
    plt.ylabel(r'$y$ (m)')    
    plt.xlim((coords_x.min(), coords_x.max()))
    plt.ylim((coords_y.min(), coords_y.max()))
    h = plt.colorbar(cm.ScalarMappable(cmap=cm.summer))
    h.set_label('$HDRF_{\\rm680 nm}$ (%s)' % (R_units), rotation=90)
    plt.savefig('1_HDRF_680.png', dpi=300)
    
    fig = plt.figure(figsize =(10.5, 9))
    ax = plt.axes()
    ax.pcolormesh(coords_x, coords_y, F_[:, :, 1], shading='auto',
                  cmap=cm.summer)
    plt.plot(FP_.Xutm, FP_.Yutm, 'or', markersize=11 )    
    for i_ in range(len(FP_)):
        plt.text(FP_.Xutm[i_] + 3*pixel_size, FP_.Yutm[i_] + 3*pixel_size,
                 ('%d' % FP_.PlotNum[i_]), color='r', fontsize=15,
                 fontweight='bold')
    plt.plot(FPhh_.XNPQwheat[0], FPhh_.YNPQwheat[0], 'sb', markersize=11)  
    plt.text(FPhh_.XNPQwheat[0] - 20*pixel_size, FPhh_.YNPQwheat[0],
                 'moni-PAM', color='b', fontsize=15, fontweight='bold')
    plt.plot(FPhh_.XNPQmaize[0], FPhh_.YNPQmaize[0], 'sb', markersize=11)  
    plt.text(FPhh_.XNPQmaize[0] + 3*pixel_size, FPhh_.YNPQmaize[0],
                 'moni-PAM', color='b', fontsize=15, fontweight='bold')
    plt.text(coords_x[0, 0] + n_cols/4,coords_y[0, 0] + 5,
        r'Field 1 (${Triticum}$ ${aestivum}$)', horizontalalignment='center',
        fontsize=12)
    plt.text(coords_x[0, 0] + 3*n_cols/4,coords_y[0, 0] + 5,
        r'Field 2 (${Zea}$ ${mays}$)', horizontalalignment='center',
        fontsize=12)
    plt.xlabel(r'$x$ (m)')
    plt.ylabel(r'$y$ (m)')    
    plt.xlim((coords_x.min(), coords_x.max()))
    plt.ylim((coords_y.min(), coords_y.max()))
    h = plt.colorbar(cm.ScalarMappable(cmap=cm.summer))
    h.set_label('$F_{\\rm760 nm}$ ($\\rm%s$)' % (F_units), rotation=90)
    plt.savefig('2_F_760.png', dpi=300)
    
    fig = plt.figure(figsize =(10.5, 9))
    ax = plt.axes()
    ax.pcolormesh(coords_x, coords_y, LST_, shading='auto',
                  cmap=cm.summer)
    plt.plot(FP_.Xutm, FP_.Yutm, 'or', markersize=11 )    
    for i_ in range(len(FP_)):
        plt.text(FP_.Xutm[i_] + 3*pixel_size, FP_.Yutm[i_] + 3*pixel_size,
                 ('%d' % FP_.PlotNum[i_]), color='r', fontsize=15,
                 fontweight='bold')
    plt.plot(FPhh_.XNPQwheat[0], FPhh_.YNPQwheat[0], 'sb', markersize=11)  
    plt.text(FPhh_.XNPQwheat[0] - 20*pixel_size, FPhh_.YNPQwheat[0],
                 'moni-PAM', color='b', fontsize=15, fontweight='bold')
    plt.plot(FPhh_.XNPQmaize[0], FPhh_.YNPQmaize[0], 'sb', markersize=11)  
    plt.text(FPhh_.XNPQmaize[0] + 3*pixel_size, FPhh_.YNPQmaize[0],
                 'moni-PAM', color='b', fontsize=15, fontweight='bold')
    plt.text(coords_x[0, 0] + n_cols/4,coords_y[0, 0] + 5,
        r'Field 1 (${Triticum}$ ${aestivum}$)', horizontalalignment='center',
        fontsize=12)
    plt.text(coords_x[0, 0] + 3*n_cols/4,coords_y[0, 0] + 5,
        r'Field 2 (${Zea}$ ${mays}$)', horizontalalignment='center',
        fontsize=12)
    plt.xlabel(r'$x$ (m)')
    plt.ylabel(r'$y$ (m)')    
    plt.xlim((coords_x.min(), coords_x.max()))
    plt.ylim((coords_y.min(), coords_y.max()))
    h = plt.colorbar(cm.ScalarMappable(cmap=cm.summer))
    h.set_label('$LST$ ($\\rm%s$)' % (F_units), rotation=90)
    plt.savefig('3_LST.png', dpi=300)

#%% 3) Brief example of how to link field and imagery data
# Generate indices for each crop
Iw = FP_['Crop']=='Wheat'
Im = FP_['Crop']!='Wheat'

# Convert row and column indices in a linear subindice
x_indices_w = FP_['Xpix'][Iw].values-1
y_indices_w = FP_['Ypix'][Iw].values-1
x_indices_m = FP_['Xpix'][Im].values-1
y_indices_m = FP_['Ypix'][Im].values-1

# # Access to the first layer of the [row, col, band] imagery matrix
plt.figure()
plt.plot(F_[x_indices_w,y_indices_w, 0],FP_['LAI'][Iw],'o',label='Wheat')
plt.plot(F_[x_indices_m,y_indices_m, 0],FP_['LAI'][Im],'s',label='Maize')
plt.xlabel('$F_{687 nm}$ ($'+F_units+'$)')
plt.ylabel('$LAI$ (m$^2$ m$^{-2}$)')
plt.grid()
if do_plots == True:
    plt.savefig('4_F687_vs_LAI.png',dpi=300,bbox_inches='tight')

# # Access to the second layer of the [row, col, band] imagery matrix
plt.figure()
plt.plot(F_[x_indices_w,y_indices_w, 1],FP_['LAI'][Iw],'o',label='Wheat')
plt.plot(F_[x_indices_m,y_indices_m, 1],FP_['LAI'][Im],'s',label='Maize')
plt.xlabel('$F_{760 nm}$ ($'+F_units+'$)')
plt.ylabel('$LAI$ (m$^2$ m$^{-2}$)')
plt.grid()
if do_plots == True:
    plt.savefig('5_F760_vs_LAI.png',dpi=300,bbox_inches='tight')

# # Access to the n-layer layer of the [row, col, band] imagery matrix
bsel_ = np.where(abs(R_wvl-680) == min(abs(R_wvl-680)))[0][0]
plt.figure()
plt.plot(R_[x_indices_w,y_indices_w, bsel_],FP_['LAI'][Iw],'o',label='Wheat')
plt.plot(R_[x_indices_m,y_indices_m, bsel_],FP_['LAI'][Im],'s',label='Maize')
plt.xlabel('$HDRF_{680 nm}$ ($'+R_units+'$)')
plt.ylabel('$LAI$ (m$^2$ m$^{-2}$)')
plt.grid()
if do_plots == True:
    plt.savefig('6_HDRF680_vs_LAI.png',dpi=300,bbox_inches='tight')

###########################################################################
# # YOUR WORK STARTS HERE!

#%% 4) Your turn. Solve the challenge!



#%% 5) Your turn. Describe your methods
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
out_methods = ori_SSCresults+'SSC_report.docx'
# out_methods = ori_SSCresults+'SSC_report.doc'
# out_methods = ori_SSCresults+'SSC_report.odt'

## 6) Your turn. Prepare your results
# # 6.1) Fill up this dictionary with your personal data
pdata_ = {}
pdata_['name'] = ''
pdata_['middle_name'] = ''
pdata_['surname'] = ''
pdata_['orcid'] = ''
pdata_['institution'] = ''
pdata_['department'] = ''
pdata_['address'] = ''
pdata_['email'] = ''
pdata_['positionexperience'] = ''
# This is the name that will be used for the compressed zip folder that 
# will be generate. Normally, it should just be the surname followed by
# the name. However, if your surname or your name included characters
# that could not be used for a folder name, write a valid version instead.
pdata_['surname_name4filename'] = pdata_['surname']+pdata_['name']

# # 6.2) Fill up this structure with your estimates and, if you did it, 
# with the estimated uncertainties of your maps. Important, these
# variables must have the same dimensions [n_rows, n_cols] than the imagery
# provided, and in the units decribed below
results = {}

# # Estimates
results['LAI_est'] = -999.*np.ones([n_rows,n_cols]) # Estimated map of leaf area index [m^2 m^-2]
results['Cab_est'] = -999.*np.ones([n_rows,n_cols]) # Estimated map of leaf chlorophyll content [ug cm^-2]
results['Vcmax25_est'] = -999.*np.ones([n_rows,n_cols]) # Estimated map of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]
results['NPQ_est'] = -999.*np.ones([n_rows,n_cols]) # Estimated map of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]

# # Stress map range between 0 for minimum stress and 1 for maximum stress (leave -999. if you did not estimate them)
results['LAI_unc'] = -999.*np.ones([n_rows,n_cols]) # Estimated uncertainties of leaf area index [m^2 m^-2]
results['Cab_unc'] = -999.*np.ones([n_rows,n_cols]) # Estimated uncertainties of leaf chlorophyll content [ug cm^-2]
results['Vcmax25_unc'] = -999.*np.ones([n_rows,n_cols]) # Estimated uncertainties of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]
results['NPQ_unc'] = -999.*np.ones([n_rows,n_cols]) # Estimated uncertainties of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]

# # Stress map range between 0 for minimum stress and 1 for maximum stress (leave -999. if you did not estimate them)
results['stress_est'] = -999.*np.ones([n_rows,n_cols]) #  % Estimated map of stress [-]

# # PERFECT!!
# Now wait for the zip file to be generated and send it by email to the
# email provided in the contact section <scalingchallenge@gmail.com>

###########################################################################
# # YOUR WORK ALMOST FINISHES HERE!

#%% 7) Spatial Scaling Challenge. Prepare standard output files
# # Check if the methods section is included the results folder
if pdata_['surname_name4filename']=='' or (not exists(out_methods)):
    warnings.warn('The results are not still ready. The zip file will not be produced')
else:
    # # Remove first the file if it exists
    if exists(out_netcdf):
        remove(out_netcdf)

    # # Create file and variables
    x = nc.Dataset(out_netcdf,'w',format='NETCDF4')
    x.createDimension('n_rows',n_rows)
    x.createDimension('n_cols',n_cols)
        
    # # Add results
    for item in results.items():
        var_name = item[0]
        var_data = item[1]
        var_nc = x.createVariable(var_name,'double',('n_rows','n_cols'))
        var_nc[:,:] = var_data[:,:]

    # # Add personal data
    for item in pdata_.items():
        x.setncattr(item[0],item[1])
    x.setncattr('Scrip_SSC', 'Python')
    # # Close file
    x.close()
       
    # # Compress the files to be sent to SENSECO Working Group 1
    zip_fname = ori_SSCresults + pdata_['surname_name4filename'] + '_4sub.zip'
    
    zipObj = ZipFile(zip_fname, 'w')
    zipObj.write(out_netcdf)
    zipObj.write(out_methods)
    zipObj.close()
    # shutil.make_archive(pdata_['surname_name4filename']+'_4subP','zip','3_SSC_results/')

    # # Final instructions
    print('Congratulations!! Your results have been stored in the following zip file:\n')
    print('\t', zip_fname,'\n\n')
    print('Send this file via email to the Spatial Scaling Challenge email adress\n')
    print('\t\t<scalingchallenge@gmail.com>\n\n')
    print('Thanks a lot for your participation!\n')
    print('Looking forward to learn from the community and share the manuscript draft any soon!\n')
    print('COST Action SENSECO Working Group 1\n')

