%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%%                     SPATIAL SCALING CHALLENGE                     %%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ------------------------------------------------------------------- % %
% %                             ORGANIZED BY                            % %
% %                                                                     % %
% %  COST ACTION CA17134 "Optical synergies for spatiotemporal SENsing  % %
% %            of Scalable ECOphysiological traits" (SENSECO)           % %
% %                        https://www.senseco.eu                       % %
% %                                                                     % %
% %                            Working Group 1                          % %
% % Closing the scaling gap: from leaf measurements to satellite images % %
% %        https://www.senseco.eu/working-groups/wg1-scaling-gap        % %
% %                                                                     % %
% %     * Dr Javier Pacheco-Labrador (vice-leader)                      % %
% %             Max Planck Institute for Biogeochemistry, Germany       % %
% %     * Dr Ma. Pilar Cendrero-Mateo (leader)                          % %
% %             University of Valencia, Spain                           % %
% %     * Dr Shari Van Wittenberghe, (vice-leader)                      % %
% %             University of Valencia, Spain                           % %
% %     * Dr Gerbrand Koren                                             % %
% %             Utrecht University, Netherlands                         % %
% %     * Prof. Zbynek Malenovský                                       % %
% %             University of Bonn, Germany                             % %
% %                                                                     % %
% % ------------------------------------------------------------------- % %
% %                           CONTACT DETAILS                           % %
% %                                                                     % %
% %   Please, contact us via <scalingchallenge@gmail.com> for any       % %
% %   question or trouble found with the code or the data provided for  % %
% %   the Spatial Scaling Challenge.                                    % %
% %   See the document SSC_description_and_instructions.pdf for further % %
% %   details about addressing questoins to the organizers.             % %
% %                                                                     % %
% % ------------------------------------------------------------------- % %
% %                      DESCRIPTION AND DISCLAIMER                     % %
% %                                                                     % %
% %   This script opens the netCDF4 and CSV files provided together     % %
% %   this code containing the datasets necessary to participate in     % %
% %   the Spatial Scaling Challenge organized by the COST Action        % %
% %   SENSECO. It also allows exporting the participant's results to    % %
% %   the only standardized format that will be accepted for submission % %
% %                                                                     % %
% %   This program is free software: you can redistribute it and/or     % %
% %   modify it under the terms of the GNU General Public License as    % % 
% %   published by the Free Software Foundation, either version 3 of    % %
% %   the License, or any later version.                                % %
% %                                                                     % %
% %   This program is distributed in the hope that it will be useful,   % %
% %   but WITHOUT ANY WARRANTY; without even the implied warranty of    % %
% %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     % %
% %   GNU General Public License for more details.                      % %
% %                                                                     % %
% %   You should have received a copy of the GNU General Public License % %
% %   along with this program. If not, see                              % %
% %   <http://www.gnu.org/licenses/>.                                   % %
% %                                                                     % %
% %   Meteorological data from Majadas de Tiétar experimental station   % %
% %   were provided by the Max Planck Institute for Biogeochemistry     % %
% %   (Germany) and Fundación CEAM (Spain)                              % %
% %                                                                     % %
% % ------------------------------------------------------------------- % %
% %                             CODE AUTHORS                            % %
% %                                                                     % %
% %     * Dr Javier Pacheco-Labrador                                    % %
% %             Max Planck Institute for Biogeochemistry, Germany       % %
% %                                                                     % %
% % ------------------------------------------------------------------- % %
% % ------------------------------------------------------------------- % %


%% 1) Initialize
close all
clear all
clc

% % Options
do_plots = true;

% % Define paths
ori_SSCdata = [pwd,'/1_SSC_data/'];
SSC_R = [ori_SSCdata, 'Airborne_HDRF.nc'];
SSC_F = [ori_SSCdata, 'Airborne_F.nc'];
SSC_LST = [ori_SSCdata, 'Airborne_LST.nc'];
SSC_field_plots = [ori_SSCdata, 'FieldData.csv'];
SSC_field_time_series = [ori_SSCdata, 'FieldData_hh.csv'];
ori_SSCresults = [pwd, '/3_SSC_results/'];
out_netcdf = [ori_SSCresults, 'maps_estimates.nc'];

%% 2) Spatial Scaling Challenge. Data import
% % Coordinates, common to all imagery.
spat_ref = ncreadatt(SSC_R, '/', 'spat_ref');
pixel_size = ncreadatt(SSC_R, '/', 'pixel_size'); %in meters
coords_x = ncread(SSC_R, 'coords_x') + pixel_size/2; %from top-right corner to center of the pixel
coords_y = ncread(SSC_R, 'coords_y') - pixel_size/2; %from top-right corner to center of the pixel
coords_z = ncread(SSC_R, 'coords_z');
[n_rows, n_cols] = size(coords_x);

% % Airborne hyperspectral reflectance imagery
R_ = ncread(SSC_R, 'var_'); %(Reflectance imagery [n_rows, n_cols, n_bands])
R_wvl = ncread(SSC_R, 'wvl_center'); %(Center wavelength)
R_fwhm = ncread(SSC_R, 'wvl_fwhm'); %(Full Width Half Maximum)
R_sza = ncread(SSC_R, 'VZA'); %View Zenith Angle
R_svaa = ncread(SSC_R, 'SVAA'); %Sun-View Azimuth Angle
R_description = ncreadatt(SSC_R, '/', 'var_description');
R_units = ncreadatt(SSC_R, '/', 'var_uds');
R_wvl_units = ncreadatt(SSC_R, '/', 'wvl_uds');
n_bands = length(R_wvl);

% % Airborne sun-induced chlorophyll fluorescence radiance imagery
F_ = ncread(SSC_F, 'var_'); %(Fluorescence radiance imagery [n_rows, n_cols, n_bands])
F_wvl = ncread(SSC_F, 'wvl_center'); %(Center wavelength)
F_fwhm = ncread(SSC_F, 'wvl_fwhm'); %(Full Width Half Maximum)
F_sza = ncread(SSC_F, 'VZA'); %View Zenith Angle
F_svaa = ncread(SSC_F, 'SVAA'); %Sun-View Azimuth Angle
F_description = ncreadatt(SSC_F, '/', 'var_description');
F_units = ncreadatt(SSC_F, '/', 'var_uds');
F_wvl_units = ncreadatt(SSC_F, '/', 'wvl_uds');

% % Airborne land surface temperature imagery
LST_ = ncread(SSC_LST, 'var_'); %(Land surface temperature imagery [n_rows, n_cols])
LST_sza = ncread(SSC_LST, 'VZA'); %View Zenith Angle
LST_svaa = ncread(SSC_LST, 'SVAA'); %Sun-View Azimuth Angle
LST_description = ncreadatt(SSC_LST, '/', 'var_description');
LST_units = ncreadatt(SSC_LST, '/', 'var_uds');

% % Field data. Spatial sampling in 1 x 1 meter plots.
FP_ = readtable(SSC_field_plots, 'delimiter', ',');

% % Field data. Meteorological data and time series of NPQ.
FPhh_ = readtable(SSC_field_time_series, 'delimiter', ',');

if do_plots;    
    % % Airborne hyperspectral reflectance imagery
    bsel_ = find(abs(R_wvl-680) == min(abs(R_wvl-680)));
    figure(1); clf; set(gcf, 'Color', 'w'); hold on; colormap('summer');
    surf(coords_x, coords_y, R_(:, :, bsel_), 'EdgeColor','none');
    p = plot3(FP_.Xutm, FP_.Yutm, ones(size(FP_.Yutm)), 'or', 'MarkerFaceColor','r');
    for i_ = 1:size(FP_, 1)
        text(FP_.Xutm(i_) + 3*pixel_size, FP_.Yutm(i_) + 3*pixel_size, 1, sprintf('%d',...
            FP_.PlotNum(i_)), 'Color', 'r', 'FontWeight', 'Bold')
    end
    phh = plot3(FPhh_.XNPQwheat, FPhh_.YNPQwheat, ones(size(FPhh_.YNPQwheat)), ...
        'sc', 'MarkerFaceColor','b');
    text(FPhh_.XNPQwheat(1) - 3*pixel_size, FPhh_.YNPQwheat(1), 1, 'moni-PAM',...
        'HorizontalAlignment','right', 'Color', 'c', 'FontWeight', 'Bold')
    plot3(FPhh_.XNPQmaize, FPhh_.YNPQmaize, ones(size(FPhh_.YNPQmaize)), ...
        'sc', 'MarkerFaceColor','b');
    text(FPhh_.XNPQmaize(1) + 3*pixel_size, FPhh_.YNPQmaize(1), 1, 'moni-PAM',...
        'HorizontalAlignment','left', 'Color', 'c', 'FontWeight', 'Bold')
    text(round(coords_x(1,1) + n_cols/4),coords_y(1,1) + 5, 1,...
        sprintf('Field 1\n(\\itTriticum aestivum\\rm)'),...
        'HorizontalAlignment', 'center');
    text(round(coords_x(1,1) + 3*n_cols/4),coords_y(1,1) + 5, 1,...
        sprintf('Field 2\n(\\itZea mays\\rm)'),...
        'HorizontalAlignment', 'center');
    xlabel('\itx\rm (m)'); ylabel('\ity\rm (m)');    
    axis equal
    axis([coords_x(1), coords_x(end), coords_y(end), coords_y(1) ])
    h = colorbar(); h.Label.String = sprintf('\\itHDRF\\rm_{680 nm} (%s)', R_units);
    saveas(figure(1),'1_HDRF_680','png')
    
    % % Airborne sun-induced clorophyll fluorescence radiance imagery
    figure(2); clf; set(gcf, 'Color', 'w'); hold on;  colormap('summer');
    surf(coords_x, coords_y, F_(:, :, 2), 'EdgeColor','none');
    plot3(FP_.Xutm, FP_.Yutm, 5*ones(size(FP_.Yutm)), 'or', 'MarkerFaceColor','r');
    for i_ = 1:size(FP_, 1)
        text(FP_.Xutm(i_) + 3*pixel_size, FP_.Yutm(i_) + 3*pixel_size, 5, sprintf('%d',...
            FP_.PlotNum(i_)), 'Color', 'r', 'FontWeight', 'Bold')
    end
    plot3(FPhh_.XNPQwheat, FPhh_.YNPQwheat, 5*ones(size(FPhh_.YNPQwheat)), ...
        'sc', 'MarkerFaceColor','b');
    text(FPhh_.XNPQwheat(1) - 3*pixel_size, FPhh_.YNPQwheat(1), 5, 'moni-PAM',...
        'HorizontalAlignment','right', 'Color', 'b', 'FontWeight', 'Bold')
    plot3(FPhh_.XNPQmaize, FPhh_.YNPQmaize, 5*ones(size(FPhh_.YNPQmaize)), ...
        'sc', 'MarkerFaceColor','b');
    text(FPhh_.XNPQmaize(1) + 3*pixel_size, FPhh_.YNPQmaize(1), 5, 'moni-PAM',...
        'HorizontalAlignment','left', 'Color', 'b', 'FontWeight', 'Bold')
    text(round(coords_x(1,1) + n_cols/4),coords_y(1,1) + 5, 5,...
        sprintf('Field 1\n(\\itTriticum aestivum\\rm)'),...
        'HorizontalAlignment', 'center');
    text(round(coords_x(1,1) + 3*n_cols/4),coords_y(1,1) + 5, 5,...
        sprintf('Field 2\n(\\itZea mays\\rm)'),...
        'HorizontalAlignment', 'center');
    xlabel('\itx\rm (m)'); ylabel('\ity\rm (m)');    
    axis equal
    axis([coords_x(1), coords_x(end), coords_y(end), coords_y(1)])    
    h = colorbar(); h.Label.String = sprintf('\\itF\\rm_{760 nm} (%s)', F_units);
    saveas(figure(2),'2_F_760','png')
    
    % % Airborne sun-induced clorophyll fluorescence radiance imagery
    figure(3); clf; set(gcf, 'Color', 'w'); hold on;  colormap('summer');
    surf(coords_x, coords_y, LST_, 'EdgeColor','none');
    p = plot3(FP_.Xutm, FP_.Yutm, 400*ones(size(FP_.Yutm)), 'or', 'MarkerFaceColor','r');
    for i_ = 1:size(FP_, 1)
        text(FP_.Xutm(i_) + 3*pixel_size, FP_.Yutm(i_) + 3*pixel_size, 400, sprintf('%d',...
            FP_.PlotNum(i_)), 'Color', 'r', 'FontWeight', 'Bold')
    end
    phh = plot3(FPhh_.XNPQwheat, FPhh_.YNPQwheat, 400*ones(size(FPhh_.YNPQwheat)), ...
        'sc', 'MarkerFaceColor','b');
    text(FPhh_.XNPQwheat(1) - 3*pixel_size, FPhh_.YNPQwheat(1), 400, 'moni-PAM',...
        'HorizontalAlignment','right', 'Color', 'b', 'FontWeight', 'Bold')
    plot3(FPhh_.XNPQmaize, FPhh_.YNPQmaize, 400*ones(size(FPhh_.YNPQmaize)), ...
        'sc', 'MarkerFaceColor','b');
    text(FPhh_.XNPQmaize(1) + 3*pixel_size, FPhh_.YNPQmaize(1), 400, 'moni-PAM',...
        'HorizontalAlignment','left', 'Color', 'b', 'FontWeight', 'Bold')
    text(round(coords_x(1,1) + n_cols/4),coords_y(1,1) + 5, 400,...
        sprintf('Field 1\n(\\itTriticum aestivum\\rm)'),...
        'HorizontalAlignment', 'center');
    text(round(coords_x(1,1) + 3*n_cols/4),coords_y(1,1) + 5, 400,...
        sprintf('Field 2\n(\\itZea mays\\rm)'),...
        'HorizontalAlignment', 'center');
    xlabel('\itx\rm (m)'); ylabel('\ity\rm (m)');    
    axis equal
    axis([coords_x(1), coords_x(end), coords_y(end), coords_y(1)])    
    h = colorbar(); h.Label.String = sprintf('\\itLST\\rm (%s)', LST_units);
    saveas(figure(3),'3_LST','png')
end

%% 3) Brief example of how to link field and imagery data
% Convert row and column indices in a linear subindice
ind_ = sub2ind([n_rows, n_cols],FP_.Ypix,FP_.Xpix);

% Generate indices for each crop
Iw = strcmp(FP_.Crop, 'Wheat');
Im = ~Iw;

% % Access to the first layer of the [row, col, band] imagery matrix
figure(4); clf; set(gcf, 'Color', 'w'); hold on; grid on;
plot(F_(ind_(Iw)), FP_.LAI(Iw), 'o')
plot(F_(ind_(Im)), FP_.LAI(Im), 's')
xlabel(sprintf('\\itF\\rm_{687 nm} (%s)', F_units))
ylabel('\itLAI\rm (m^2 m^{-2} )')
if do_plots
    saveas(figure(3),'4_F687_vs_LAI','png')
end

% % Access to the second layer of the [row, col, band] imagery matrix
figure(5); clf; set(gcf, 'Color', 'w'); hold on; grid on;
plot(F_(n_rows*n_cols + ind_(Iw)), FP_.LAI(Iw), 'o')
plot(F_(n_rows*n_cols + ind_(Im)), FP_.LAI(Im), 's')
xlabel(sprintf('\\itF\\rm_{760 nm} (%s)', F_units))
ylabel('\itLAI\rm (m^2 m^{-2} )');
if do_plots
    saveas(figure(5),'5_F760_vs_LAI','png')
end

% % Access to the n-layer layer of the [row, col, band] imagery matrix
bsel_ = find(abs(R_wvl-680) == min(abs(R_wvl-680)));
figure(6); clf; set(gcf, 'Color', 'w'); hold on; grid on;
plot(R_((bsel_-1).*n_rows*n_cols + ind_(Iw)), FP_.LAI(Iw), 'o')
plot(R_((bsel_-1).*n_rows*n_cols + ind_(Im)), FP_.LAI(Im), 's')
xlabel(sprintf('\\itHDRF\\rm_{680 nm} (%s)', R_units))
ylabel('\itLAI\rm (m^2 m^{-2} )');
if do_plots
    saveas(figure(6),'6_HDRF680_vs_LAI','png')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % YOUR WORK STARTS HERE!

%% 4) Your turn. Solve the challenge!



%% 5) Your turn. Describe your methods
% Before submitting your results, you need to document the methods you
% used in a standardized way that will help to summarize the contribution
% of all the participants. Use any of the MS Word (.doc or .docx) or 
% OpenOffice Writer (.odt) templates: /2_SSC_templates/SSC_report.xxx
% First, provide your diagnosis of the actual vegetation status. Then,
% describe the methods you used to estimate each of the variables (and
% uncertainties if you did). Try to be concise and clear. The descriptions
% could be included in the supplementary material of the joint manuscript,
% therefore, take care of English grammar and style.
% Include references if necessary (Author et al., year) and the full
% reference at the end of each section.

% % Here, select the filename matching the extension of your methods' file
out_methods = [ori_SSCresults,'SSC_report.docx'];
% out_methods = [ori_SSCresults,'SSC_report.doc'];
% out_methods = [ori_SSCresults,'SSC_report.odt'];

%% 6) Your turn. Prepare your results
% % 6.1) Fill up this structure with your personal data
pdata_ = struct();
pdata_.name = '';
pdata_.middle_name = '';
pdata_.surname = '';
pdata_.orcid = '';
pdata_.institution = '';
pdata_.department = '';
pdata_.address = '';
pdata_.email = '';
pdata_.positionexperience = '';
% This is the name that will be used for the compressed zip folder that 
% will be generate. Normally, it should just be the surname followed by
% the name. However, if your surname or your name included characters
% that could not be used for a folder name, write a valid version instead.
pdata_.surname_name4filename = [pdata_.surname,pdata_.name];

% % 6.2) Fill up this structure with your estimates and, if you did it, 
% with the estimated uncertainties of your maps. Important, these
% variables must have the same dimensions [n_rows, n_cols] than the imagery
% provided, and in the units decribed below
results = struct();

% % Estimates
results.LAI_est = -999.*ones(n_rows, n_cols); % Estimated map of leaf area index [m^2 m^-2]
results.Cab_est = -999.*ones(n_rows, n_cols); % Estimated map of leaf chlorophyll content [ug cm^-2]
results.Vcmax25_est = -999.*ones(n_rows, n_cols); % Estimated map of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]
results.NPQ_est = -999.*ones(n_rows, n_cols); % Estimated map of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]

% % Uncertainties (leave -999. if you did not estimate them)
results.LAI_unc = -999.*ones(n_rows, n_cols); % Estimated uncertainties of leaf area index [m^2 m^-2]
results.Cab_unc = -999.*ones(n_rows, n_cols); % Estimated uncertainties of leaf chlorophyll content [ug cm^-2]
results.Vcmax25_unc = -999.*ones(n_rows, n_cols); % Estimated uncertainties of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]
results.NPQ_unc = -999.*ones(n_rows, n_cols); % Estimated uncertainties of maximum carboxylation rate at 25 ºC [umol cm^-2 s^-1]

% % Stress map range between 0 for minimum stress and 1 for maximum stress (leave -999. if you did not estimate them) %%
results.stress_est = -999.*ones(n_rows, n_cols); % Estimated map of stress [-] %%

% % PERFECT!!
% Now wait for the zip file to be generated and send it by email to the
% email provided in the contact section <scalingchallenge@gmail.com>


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % YOUR WORK ALMOST FINISHES HERE!

%% 7) Spatial Scaling Challenge. Prepare standar output files
% % Check if the methods section is included the results folder
if strcmp(pdata_.surname_name4filename, '') || ~exist(out_methods, 'file')
    warning('The results are not still ready. The zip file will not be produced')
else
% % Remove first the file if it exists
    if exist(out_netcdf, 'file')
        delete(out_netcdf)
    end
    
% % Add results
    fld  = fieldnames(results);
    for i_ = 1:length(fld)        
        nccreate(out_netcdf, fld{i_}, 'Dimensions',...
            {'rows', size(results.(fld{i_}), 1),...
            'cols', size(results.(fld{i_}), 1)}, 'Datatype', 'double',...
            'Format', 'netcdf4', 'FillValue', -99.)
        ncwrite(out_netcdf, fld{i_},  results.(fld{i_}));
    end
    
% % Add personal data
    fld = fieldnames(pdata_);
    for i_ = 1:length(fld)
        ncwriteatt(out_netcdf, '/', fld{i_}, pdata_.(fld{i_}))
    end
    ncwriteatt(out_netcdf, '/', 'Scrip_SSC', 'Matlab')
    
% % Compress the files to be sent to SENSECO Working Group 1
    zip_fname = [ori_SSCresults, pdata_.surname_name4filename, '_4sub.zip']; 
    zip(zip_fname, {out_netcdf, out_methods})
    
% % Final instructions
    fprintf('Congratulations!! Your results have been stored in the following zip file:\n');   
    fprintf('\t%s\n\n', zip_fname);
    fprintf('Send this file via email to the Spatial Scaling Challenge email adress\n');
    fprintf('\t\t<scalingchallenge@gmail.com>\n\n');
    fprintf('Thanks a lot for your participation!\n');
    fprintf('Looking forward to learn from the community and share the manuscript draft any soon!\n');
    fprintf('COST Action SENSECO Working Group 1\n');
end