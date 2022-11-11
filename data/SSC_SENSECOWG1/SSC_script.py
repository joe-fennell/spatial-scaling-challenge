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
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import cm
# # import netCDF4 as nc
# import pandas as pd
from os.path import exists
from os import remove
import os
import warnings
from zipfile import ZipFile

# # # Options
# do_plots = True

# # # Define paths
# ori_SSCdata = '1_SSC_data//'
# SSC_R = ori_SSCdata + 'Airborne_HDRF.nc'
# SSC_F = ori_SSCdata + 'Airborne_F.nc'
# SSC_LST = ori_SSCdata + 'Airborne_LST.nc'
# SSC_field_plots = ori_SSCdata + 'FieldData.csv'
# SSC_field_time_series = ori_SSCdata + 'FieldData_hh.csv'
ori_SSCresults = '3_SSC_results'
out_netcdf = os.path.join(ori_SSCresults, 'maps_estimates.nc')

# #%% 2) Spatial Scaling Challenge. Data import
# # # Coordinates, common to all imagery. 
# h = nc.Dataset(SSC_R)
# spat_ref = h.getncattr('spat_ref')
# pixel_size = h.getncattr('pixel_size') #in meters
# coords_x =  h.variables['coords_x'][:] + pixel_size/2 #from top-right corner to center of the pixel
# coords_y =  h.variables['coords_y'][:] - pixel_size/2 #from top-right corner to center of the pixel
# coords_z =  h.variables['coords_z'][:]
# [n_rows, n_cols] = np.shape(coords_x)
# h.close()

# # # Airborne hyperspectral reflectance imagery
# h = nc.Dataset(SSC_R)
# R_ = np.transpose(h.variables['var_'][:],(1, 2, 0)) #(Reflectance imagery [n_rows, n_cols, n_bands])
# R_wvl = h.variables['wvl_center'][:] #(Center wavelength)
# R_fwhm = h.variables['wvl_fwhm'][:] #(Full Width Half Maximum)
# R_sza = h.variables['VZA'] #(View Zenith Angle)
# R_svaa = h.variables['SVAA'] #(Sun-View Azimuth Angle)
# R_description = h.getncattr('var_description')
# R_units = h.getncattr('var_uds')
# R_wvl_units = h.getncattr('wvl_uds')
# n_bands = len(R_wvl)
# h.close()

# # # Airborne sun-induced clorophyll fluorescence radiance imagery
# g = nc.Dataset(SSC_F)
# F_ = np.transpose(g.variables['var_'][:],(1, 2, 0)) #(Fluorescence radiance imagery [n_rows, n_cols, n_bands])
# F_wvl = g.variables['wvl_center'][:] #(Center wavelength)
# F_fwhm = g.variables['wvl_fwhm'][:] #(Full Width Half Maximum)
# F_sza = g.variables['VZA'] #(View Zenith Angle)
# F_svaa = g.variables['SVAA'] #(Sun-View Azimuth Angle)
# F_description = g.getncattr('var_description')
# F_units = g.getncattr('var_uds')
# F_wvl_units = g.getncattr('wvl_uds')
# g.close()

# # # Airborne land surface temperature imagery
# q = nc.Dataset(SSC_LST)
# LST_ = np.transpose(q.variables['var_'][:],(1, 2, 0))[:, :, 0] #(Land surface temperature imagery, [n_rows, n_cols])
# LST_sza = q.variables['VZA'] #(View Zenith Angle)
# LST_svaa = q.variables['SVAA'] #(Sun-View Azimuth Angle)
# LST_description = q.getncattr('var_description')
# LST_units = q.getncattr('var_uds')
# q.close()

# # # Field data. Spatial sampling in 1 x 1 meter plots.
# FP_ = pd.read_csv(SSC_field_plots, delimiter=',')

# # # Field data. Meteorological data and time series of NPQ.
# FPhh_ = pd.read_csv(SSC_field_time_series, delimiter=',')

# if do_plots is True:    
#     # # Airborne hyperspectral reflectance imagery
#     bsel_ = np.where(abs(R_wvl-680) == min(abs(R_wvl-680)))[0][0]
#     fig = plt.figure(figsize =(10.5, 9))
#     ax = plt.axes()
#     ax.pcolormesh(coords_x, coords_y, R_[:, :, bsel_], shading='auto',
#                   cmap=cm.summer)
#     plt.plot(FP_.Xutm, FP_.Yutm, 'or', markersize=11 )    
#     for i_ in range(len(FP_)):
#         plt.text(FP_.Xutm[i_] + 3*pixel_size, FP_.Yutm[i_] + 3*pixel_size,
#                  ('%d' % FP_.PlotNum[i_]), color='r', fontsize=15,
#                  fontweight='bold')
#     plt.plot(FPhh_.XNPQwheat[0], FPhh_.YNPQwheat[0], 'sb', markersize=11)  
#     plt.text(FPhh_.XNPQwheat[0] - 20*pixel_size, FPhh_.YNPQwheat[0],
#                  'moni-PAM', color='b', fontsize=15, fontweight='bold')
#     plt.plot(FPhh_.XNPQmaize[0], FPhh_.YNPQmaize[0], 'sb', markersize=11)  
#     plt.text(FPhh_.XNPQmaize[0] + 3*pixel_size, FPhh_.YNPQmaize[0],
#                  'moni-PAM', color='b', fontsize=15, fontweight='bold')
#     plt.text(coords_x[0, 0] + n_cols/4,coords_y[0, 0] + 5,
#         r'Field 1 (${Triticum}$ ${aestivum}$)', horizontalalignment='center',
#         fontsize=12)
#     plt.text(coords_x[0, 0] + 3*n_cols/4,coords_y[0, 0] + 5,
#         r'Field 2 (${Zea}$ ${mays}$)', horizontalalignment='center',
#         fontsize=12)
#     plt.xlabel(r'$x$ (m)')
#     plt.ylabel(r'$y$ (m)')    
#     plt.xlim((coords_x.min(), coords_x.max()))
#     plt.ylim((coords_y.min(), coords_y.max()))
#     h = plt.colorbar(cm.ScalarMappable(cmap=cm.summer))
#     h.set_label('$HDRF_{\\rm680 nm}$ (%s)' % (R_units), rotation=90)
#     plt.savefig('1_HDRF_680.png', dpi=300)
    
#     fig = plt.figure(figsize =(10.5, 9))
#     ax = plt.axes()
#     ax.pcolormesh(coords_x, coords_y, F_[:, :, 1], shading='auto',
#                   cmap=cm.summer)
#     plt.plot(FP_.Xutm, FP_.Yutm, 'or', markersize=11 )    
#     for i_ in range(len(FP_)):
#         plt.text(FP_.Xutm[i_] + 3*pixel_size, FP_.Yutm[i_] + 3*pixel_size,
#                  ('%d' % FP_.PlotNum[i_]), color='r', fontsize=15,
#                  fontweight='bold')
#     plt.plot(FPhh_.XNPQwheat[0], FPhh_.YNPQwheat[0], 'sb', markersize=11)  
#     plt.text(FPhh_.XNPQwheat[0] - 20*pixel_size, FPhh_.YNPQwheat[0],
#                  'moni-PAM', color='b', fontsize=15, fontweight='bold')
#     plt.plot(FPhh_.XNPQmaize[0], FPhh_.YNPQmaize[0], 'sb', markersize=11)  
#     plt.text(FPhh_.XNPQmaize[0] + 3*pixel_size, FPhh_.YNPQmaize[0],
#                  'moni-PAM', color='b', fontsize=15, fontweight='bold')
#     plt.text(coords_x[0, 0] + n_cols/4,coords_y[0, 0] + 5,
#         r'Field 1 (${Triticum}$ ${aestivum}$)', horizontalalignment='center',
#         fontsize=12)
#     plt.text(coords_x[0, 0] + 3*n_cols/4,coords_y[0, 0] + 5,
#         r'Field 2 (${Zea}$ ${mays}$)', horizontalalignment='center',
#         fontsize=12)
#     plt.xlabel(r'$x$ (m)')
#     plt.ylabel(r'$y$ (m)')    
#     plt.xlim((coords_x.min(), coords_x.max()))
#     plt.ylim((coords_y.min(), coords_y.max()))
#     h = plt.colorbar(cm.ScalarMappable(cmap=cm.summer))
#     h.set_label('$F_{\\rm760 nm}$ ($\\rm%s$)' % (F_units), rotation=90)
#     plt.savefig('2_F_760.png', dpi=300)
    
#     fig = plt.figure(figsize =(10.5, 9))
#     ax = plt.axes()
#     ax.pcolormesh(coords_x, coords_y, LST_, shading='auto',
#                   cmap=cm.summer)
#     plt.plot(FP_.Xutm, FP_.Yutm, 'or', markersize=11 )    
#     for i_ in range(len(FP_)):
#         plt.text(FP_.Xutm[i_] + 3*pixel_size, FP_.Yutm[i_] + 3*pixel_size,
#                  ('%d' % FP_.PlotNum[i_]), color='r', fontsize=15,
#                  fontweight='bold')
#     plt.plot(FPhh_.XNPQwheat[0], FPhh_.YNPQwheat[0], 'sb', markersize=11)  
#     plt.text(FPhh_.XNPQwheat[0] - 20*pixel_size, FPhh_.YNPQwheat[0],
#                  'moni-PAM', color='b', fontsize=15, fontweight='bold')
#     plt.plot(FPhh_.XNPQmaize[0], FPhh_.YNPQmaize[0], 'sb', markersize=11)  
#     plt.text(FPhh_.XNPQmaize[0] + 3*pixel_size, FPhh_.YNPQmaize[0],
#                  'moni-PAM', color='b', fontsize=15, fontweight='bold')
#     plt.text(coords_x[0, 0] + n_cols/4,coords_y[0, 0] + 5,
#         r'Field 1 (${Triticum}$ ${aestivum}$)', horizontalalignment='center',
#         fontsize=12)
#     plt.text(coords_x[0, 0] + 3*n_cols/4,coords_y[0, 0] + 5,
#         r'Field 2 (${Zea}$ ${mays}$)', horizontalalignment='center',
#         fontsize=12)
#     plt.xlabel(r'$x$ (m)')
#     plt.ylabel(r'$y$ (m)')    
#     plt.xlim((coords_x.min(), coords_x.max()))
#     plt.ylim((coords_y.min(), coords_y.max()))
#     h = plt.colorbar(cm.ScalarMappable(cmap=cm.summer))
#     h.set_label('$LST$ ($\\rm%s$)' % (F_units), rotation=90)
#     plt.savefig('3_LST.png', dpi=300)

# #%% 3) Brief example of how to link field and imagery data
# # Generate indices for each crop
# Iw = FP_['Crop']=='Wheat'
# Im = FP_['Crop']!='Wheat'

# # Convert row and column indices in a linear subindice
# x_indices_w = FP_['Xpix'][Iw].values-1
# y_indices_w = FP_['Ypix'][Iw].values-1
# x_indices_m = FP_['Xpix'][Im].values-1
# y_indices_m = FP_['Ypix'][Im].values-1

# # # Access to the first layer of the [row, col, band] imagery matrix
# plt.figure()
# plt.plot(F_[x_indices_w,y_indices_w, 0],FP_['LAI'][Iw],'o',label='Wheat')
# plt.plot(F_[x_indices_m,y_indices_m, 0],FP_['LAI'][Im],'s',label='Maize')
# plt.xlabel('$F_{687 nm}$ ($'+F_units+'$)')
# plt.ylabel('$LAI$ (m$^2$ m$^{-2}$)')
# plt.grid()
# if do_plots == True:
#     plt.savefig('4_F687_vs_LAI.png',dpi=300,bbox_inches='tight')

# # # Access to the second layer of the [row, col, band] imagery matrix
# plt.figure()
# plt.plot(F_[x_indices_w,y_indices_w, 1],FP_['LAI'][Iw],'o',label='Wheat')
# plt.plot(F_[x_indices_m,y_indices_m, 1],FP_['LAI'][Im],'s',label='Maize')
# plt.xlabel('$F_{760 nm}$ ($'+F_units+'$)')
# plt.ylabel('$LAI$ (m$^2$ m$^{-2}$)')
# plt.grid()
# if do_plots == True:
#     plt.savefig('5_F760_vs_LAI.png',dpi=300,bbox_inches='tight')

# # # Access to the n-layer layer of the [row, col, band] imagery matrix
# bsel_ = np.where(abs(R_wvl-680) == min(abs(R_wvl-680)))[0][0]
# plt.figure()
# plt.plot(R_[x_indices_w,y_indices_w, bsel_],FP_['LAI'][Iw],'o',label='Wheat')
# plt.plot(R_[x_indices_m,y_indices_m, bsel_],FP_['LAI'][Im],'s',label='Maize')
# plt.xlabel('$HDRF_{680 nm}$ ($'+R_units+'$)')
# plt.ylabel('$LAI$ (m$^2$ m$^{-2}$)')
# plt.grid()
# if do_plots == True:
#     plt.savefig('6_HDRF680_vs_LAI.png',dpi=300,bbox_inches='tight')

###########################################################################
# # YOUR WORK STARTS HERE!

#%% 4) Your turn. Solve the challenge!
import xarray
import pandas as pd
import numpy as np

from sklearn.base import BaseEstimator
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.preprocessing import MinMaxScaler, RobustScaler
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn.gaussian_process.kernels import WhiteKernel, ConstantKernel, Matern

def stack(x, geometry=False):
    _stacked = x.stack(pix=['row','col'])
    if geometry:
        return np.vstack([_stacked.row.values, _stacked.col.values]).T
    return np.atleast_2d(_stacked.values)

def preprocess_hsi(x, data_var='var_',
                   geometry_vars=[]):
    xnew=[]
    xnew.append(stack(x[data_var]))
    for v in geometry_vars:
        xnew.append(stack(x[v]))
    return np.concatenate(xnew).T

def to_xarray(x, shape=(100,100,5)):
    return xarray.DataArray(x.reshape(shape).T,
                         dims=['dim0','col','row']).transpose('dim0', 'row', 'col')

def get_pix_n(xs, ys):
    lookup = np.arange(10000).reshape((100,100))
    return lookup[ys, xs]

def get_window_n(xs, ys, k=[-1,0,1]):
    """same as get_pix_n but retrieves a windowed region around point
    """
    out = []
    for i in k:
        for j in k:
            out.append(get_pix_n(xs+i, ys+j))
    return np.concatenate(out)

def to_result_array(x):
    # transpose
    x = x.transpose('row', 'col', ...)
    # fillna with -999
    x = x.fillna(-999)
    try:
        return x.values[:,:, 0]
    except IndexError:
        return x.values

class RobustScaler1D(RobustScaler):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def fit(self, X, y=None):
        super().fit(X.reshape((-1,1)), y)
    
    def transform(self, X):
        X_new = super().transform(X.reshape((-1,1)))
        return X_new.ravel()
    
    def fit_transform(self, X, y=None):
        _X = X.reshape((-1,1))
        X_new = super().fit(_X).transform(_X)
        return X_new.ravel()
    
    def inverse_transform(self, X):
        try:
            _X = X.reshape((-1,1))
            return super().inverse_transform(_X).ravel()
        except AttributeError:
            mu = super().inverse_transform(X[0].reshape((-1,1))).ravel()
            lwr = super().inverse_transform((X[0] - X[1]).reshape((-1,1))).ravel()
            SE = mu - lwr
            return mu, SE

# Open all datasets using pandas and xarray
aerial_F = xarray.open_dataset('1_SSC_data/Airborne_F.nc')
aerial_HSI = xarray.open_dataset('1_SSC_data/Airborne_HDRF.nc')
aerial_LST = xarray.open_dataset('1_SSC_data/Airborne_LST.nc')
field_data = pd.read_csv('1_SSC_data/FieldData.csv')
field_data2 = pd.read_csv('1_SSC_data/FieldData_hh.csv')       

# Apply median filtering in spectral domain
aerial_HSI_filtered = aerial_HSI['var_'].rolling({'bands':3}, center=True).median().bfill('bands').ffill('bands')

# generate a 2D dataset
N_COMPONENTS = 5

# Do a PCA decomposition
spectral_decomposer = Pipeline(
    [('rescale', RobustScaler(unit_variance=True)),
     ('pca', PCA(n_components=N_COMPONENTS))])

PCs = to_xarray(spectral_decomposer.fit_transform(
    stack(aerial_HSI_filtered).T), (100, 100, -1))

# use a threshold in the first PC to define field boundary region
field_boundary = PCs.isel(dim0=0) >= 0
field_boundary.plot()

# convert to a standard input format
gpr_input = preprocess_hsi(xarray.Dataset({'var_':PCs}))

# get the PC features and labels from field data

X = gpr_input[get_window_n(field_data['Xpix']-1, field_data['Ypix']-1, [-1,0,1]),...]

###### LAI #######
lai_scaler = RobustScaler1D()

y = lai_scaler.fit_transform(np.tile(field_data['LAI'].values, 9))

# Specify Gaussian Kernel
kernel = Matern(length_scale=np.ones(N_COMPONENTS), nu=5/2, length_scale_bounds=[.1, 1e8]) + \
WhiteKernel() 


GPR_LA = Pipeline([
     ('rescale_pca', MinMaxScaler()),
     ('regressor', GaussianProcessRegressor(kernel=kernel))]
)

# Fit the GPR
GPR_LA.fit(X, y)

# Generate predictions
_prediction_LA, _SE_LA = lai_scaler.inverse_transform(GPR_LA.predict(gpr_input, return_std=True))
prediction_LA = field_boundary.copy()
prediction_LA.values = _prediction_LA.reshape((100,100))
prediction_LA = prediction_LA.where(~field_boundary)

prediction_LA_SE = field_boundary.copy().where(~field_boundary)
prediction_LA_SE.values = _SE_LA.reshape((100,100))
prediction_LA_SE = prediction_LA_SE.where(~field_boundary)

###### Cab #######
cab_scaler = RobustScaler1D() # scale the response var for numerical stability
y2 = cab_scaler.fit_transform(np.tile(field_data['Cab'].values, 9))

# Specify Gaussian Kernel
kernel = Matern(length_scale=np.ones(N_COMPONENTS), nu=5/2) + \
WhiteKernel()

GPR_cab = Pipeline([
    #('rescale', RobustScaler()),
    #('PCA', PCA(n_components=N_COMPONENTS)),
     ('rescale_pca', MinMaxScaler()),
     ('regressor', GaussianProcessRegressor(kernel=kernel))]
)


# Fit the GPR
GPR_cab.fit(X, y2)

# Predict Cab
_prediction_cab, _SE_cab = cab_scaler.inverse_transform(GPR_cab.predict(gpr_input, return_std=True))
prediction_cab = field_boundary.copy()
prediction_cab.values = _prediction_cab.reshape((100,100))
prediction_cab = prediction_cab.where(~field_boundary)

prediction_cab_SE = field_boundary.copy().where(~field_boundary)
prediction_cab_SE.values = _SE_cab.reshape((100,100))
prediction_cab_SE = prediction_cab_SE.where(~field_boundary)

###### VCMax #######
vcmax_scaler = RobustScaler1D()

y3 = np.tile(vcmax_scaler.fit_transform(field_data['Vcmax25'].values),9)

# Specify Gaussian Kernel
kernel = Matern(length_scale=np.ones(X.shape[1]), nu=5/2) + \
WhiteKernel() 


GPR_vcmax = Pipeline(
    [('rescale_pca', RobustScaler()),
     ('regressor', GaussianProcessRegressor(kernel=kernel))]
)

# Fit the GPR
GPR_vcmax.fit(X, y3)

# predict
_prediction_vcm, _SE_vcm = vcmax_scaler.inverse_transform(GPR_vcmax.predict(gpr_input, return_std=True))
prediction_vcm = field_boundary.copy()
prediction_vcm.values = _prediction_vcm.reshape((100,100))
prediction_vcm = prediction_vcm.where(~field_boundary)

prediction_vcm_SE = field_boundary.copy().where(~field_boundary)
prediction_vcm_SE.values = _SE_vcm.reshape((100,100))
prediction_vcm_SE = prediction_vcm_SE.where(~field_boundary)

###### NPQ #######
npq_scaler = RobustScaler1D()

gpr_input2 = np.vstack([
    # stack(prediction_cab.fillna(0)),
    stack(prediction_LA.fillna(0)),
    stack(aerial_F['var_'].isel(bands=0)).ravel(),
    stack(aerial_F['var_'].isel(bands=1)).ravel()
]).T

# larger window for this due to lower GSD of SIF imager
X2 = gpr_input2[get_window_n(field_data['Xpix']-1, field_data['Ypix']-1, [-2,-1,0,1,2]),...]
y4 = np.tile(npq_scaler.fit_transform(field_data['NPQ1'].values), 25)

# Specify Gaussian Kernel
kernel = Matern(length_scale=np.ones(X2.shape[1]), nu=5/2) + \
WhiteKernel() 


GPR_npq = Pipeline(
    [('rescale_pca', RobustScaler()),
     ('regressor', GaussianProcessRegressor(kernel=kernel))]
)

# Fit the GPR
GPR_npq.fit(X2, y4)

# predict
_prediction_npq, _SE_npq = npq_scaler.inverse_transform(GPR_npq.predict(gpr_input2, return_std=True))
prediction_npq = field_boundary.copy()
prediction_npq.values = _prediction_npq.reshape((100,100))
prediction_npq = prediction_npq.where(~field_boundary)

prediction_npq_SE = field_boundary.copy().where(~field_boundary)
prediction_npq_SE.values = _SE_npq.reshape((100,100))
prediction_npq_SE = prediction_npq_SE.where(~field_boundary)
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
pdata_['name'] = 'Joseph'
pdata_['middle_name'] = 'T.'
pdata_['surname'] = 'Fennell'
pdata_['orcid'] = '0000-0001-6874-6667'
pdata_['institution'] = 'The Open University'
pdata_['department'] = 'School of Environment, Earth and Ecosystem Sciences'
pdata_['address'] = 'The Open University, Walton Hall, Milton Keynes, United Kingdom'
pdata_['email'] = 'joseph.fennell@open.ac.uk'
pdata_['positionexperience'] = 'Research Fellow'
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
n_rows = 100
n_cols = 100
results['LAI_est'] = to_result_array(prediction_LA) # Estimated map of leaf area index [m^2 m^-2]
results['Cab_est'] = to_result_array(prediction_cab) # Estimated map of leaf chlorophyll content [ug cm^-2]
results['Vcmax25_est'] = to_result_array(prediction_vcm) # Estimated map of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]
results['NPQ_est'] = to_result_array(prediction_npq) # Estimated map of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]

# # Stress map range between 0 for minimum stress and 1 for maximum stress (leave -999. if you did not estimate them)
results['LAI_unc'] = to_result_array(prediction_LA_SE) # Estimated uncertainties of leaf area index [m^2 m^-2]
results['Cab_unc'] = to_result_array(prediction_cab_SE) # Estimated uncertainties of leaf chlorophyll content [ug cm^-2]
results['Vcmax25_unc'] = to_result_array(prediction_vcm_SE) # Estimated uncertainties of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]
results['NPQ_unc'] = to_result_array(prediction_npq_SE) # Estimated uncertainties of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]

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

