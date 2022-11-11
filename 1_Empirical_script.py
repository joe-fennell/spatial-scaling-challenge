import xarray
import pandas as pd
import matplotlib.pyplot as plt
import seaborn
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
aerial_F = xarray.open_dataset('data/SSC_SENSECOWG1/1_SSC_data/Airborne_F.nc')
aerial_HSI = xarray.open_dataset('data/SSC_SENSECOWG1/1_SSC_data/Airborne_HDRF.nc')
aerial_LST = xarray.open_dataset('data/SSC_SENSECOWG1/1_SSC_data/Airborne_LST.nc')
field_data = pd.read_csv('data/SSC_SENSECOWG1/1_SSC_data/FieldData.csv')
field_data2 = pd.read_csv('data/SSC_SENSECOWG1/1_SSC_data/FieldData_hh.csv')       

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

