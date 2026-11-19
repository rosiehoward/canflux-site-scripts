% CanFlux_CFI_Figures.m
%
% Largely based on createCAFluxSiteMap.m
%
% Rosie Howard
% 18 November 2026

% Section list:
%   0. Extract data from relevant spreadsheet
%   1. Find active vs. historical vs. known
%   2a. Network (Ameriflux)
%   2b. BASE data (on Ameriflux)
%   2c. Processing (FLUXNET)


%% 0. Extract data from TSV file
clear;
fpath = '/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux';
% fname = 'AmeriFlux-sites-CanadianORIGDOWNLOAD_19August2025.tsv';    % Ameriflux data; original download
% fname_3Oct2025 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet3Oct2025_withLatLons.tsv';   % only sites with LatLons
% fname_7Oct2025 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet7Oct2025_edited.tsv';     % Ameriflux + self-reported
% fname_21Oct2025 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet21Oct2025.tsv';     % Ameriflux + self-reported
% fname_28Nov2025 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet28Nov2025.tsv';     % Ameriflux + self-reported
% fname_6Feb2026 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet6Feb2026.tsv';     % Ameriflux + self-reported
% fname_10Feb2026 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet10Feb2026_NotIncludingProposedSites.tsv';   % Ameriflux + self-reported (not including proposed sites) -> presented on Can-Peat 2026 poster
% fname_20Feb2026 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet20Feb2026_NotIncludingProposedSites.tsv';   % Ameriflux + self-reported (not including proposed sites) -> for CFI proposal

fname_18Nov2026 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet18Nov2026_NotIncludingProposedSites.tsv';   % Ameriflux + self-reported (not including proposed sites) -> for CFI proposal
w = warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
t = readtable(fullfile(fpath,fname_18Nov2026), "FileType","text",'Delimiter', '\t');                    % does not include proposed towers

fname_ORIG_17Sept2025 = 'Canadian-Flux-Sites_17Sept2025_ORIG_forComparison.tsv';
v = readtable(fullfile(fpath,fname_ORIG_17Sept2025), "FileType","text",'Delimiter', '\t');

fname_18Nov2026_incProposedSites = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet18Nov2026_IncludingProposedSites.tsv';
x = readtable(fullfile(fpath,fname_18Nov2026_incProposedSites), "FileType","text",'Delimiter', '\t');   % includes proposed towers

% fname_Baldocchi = 'CanadianSites_Baldocchi2025.csv';
% t = readtable(fullfile(fpath,fname_Baldocchi),'Delimiter', ',');

addpath(genpath('/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux'));


%% 1. Find active vs. historical vs. unknown 

activeSitesCol = t.ActiveSite;

activeUnknown = isnan(activeSitesCol);
n_activeUnknown = sum(activeUnknown) - 2;    % subtract 2 for the subheaders in column 1 (siteID)
ind_activeUnknown = find(activeUnknown);     % includes indices for rows with subheaders in column 1 (siteID)
ind = find(strcmp(t.AmerifluxSiteID,'Sites likely not registered on Ameriflux'));   % identify index dividing Ameriflux and non-Ameriflux sites in spreadsheet
ind_remove = find(ind_activeUnknown == ind);
ind_activeUnknown(ind_remove) = [];
ind_activeUnknown(1) = [];

ind_activeSites = find(activeSitesCol > 0);
n_activeSites = length(ind_activeSites);

ind_inactiveSites = find(activeSitesCol == 0);
n_inactiveSites = length(ind_inactiveSites);

ind_nocoords = find(isnan(t.Latitude_degrees_));
% ind_nocoords = find(isnan(t.Latitude_degrees_));

% remove site listings if lat/lon doesn't exist (for all categories)
common_values = intersect(ind_activeSites,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_activeSites ~= common_values(i));
    ind_activeSites = ind_activeSites(mask);
end
n_activeSites_WithCoords = length(ind_activeSites);

% remove site listings if lat/lon doesn't exist (for all categories)
common_values = intersect(ind_inactiveSites,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_inactiveSites ~= common_values(i));
    ind_inactiveSites = ind_inactiveSites(mask);
end
n_inactiveSites_WithCoords = length(ind_inactiveSites);

% remove site listings if lat/lon doesn't exist (for all categories)
common_values = intersect(ind_activeUnknown,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_activeUnknown ~= common_values(i));
    ind_activeUnknown = ind_activeUnknown(mask);
end
n_activeUnknown_WithCoords = length(ind_activeUnknown);


%% 2. Network, processing, and data (registered on Ameriflux/FLUXNET-processed/has data on Ameriflux)

% 2a. Ameriflux
[a,b] = size(t);
ind_divide = find(strcmp(t.AmerifluxSiteID,'Sites likely not registered on Ameriflux'));   % identify index dividing Ameriflux and non-Ameriflux sites in spreadsheet

ind_Ameriflux = 2:ind_divide-1;    % start at second index, until one before the dividing subheader
ind_nonAmeriflux = ind_divide+1:a;

n_Ameriflux = length(ind_Ameriflux);
n_nonAmeriflux = length(ind_nonAmeriflux);
% check
n_totalSites = n_Ameriflux + n_nonAmeriflux;

n_Ameriflux_WithCoords = length(find(~isnan(t.Latitude_degrees_(1:ind_divide-1))));
n_nonAmeriflux_WithCoords = length(find(~isnan(t.Latitude_degrees_(ind_divide+1:a))));
% check
% n_totalSites_WithCoords = n_Ameriflux_WithCoords + n_nonAmeriflux_WithCoords;

%% 2b. BASE data on Ameriflux

% these should all have coords since they will be on Ameriflux by definition
% sites with data
ind_Ameriflux_BASE = find(~isnan(t.AmeriFluxBASEDataStart));
n_Ameriflux_BASE = length(ind_Ameriflux_BASE);

ind_Ameriflux_noBASE = find(isnan(t.AmeriFluxBASEDataStart(ind_Ameriflux)));
n_Ameriflux_noBASE = length(ind_Ameriflux_noBASE);
% check
% n_BASE_noBASE = n_Ameriflux_BASE + n_Ameriflux_noBASE;

%% 2c. FLUXNET

% these should all have coords since they will be on Ameriflux by
% definition, and also all be a subset of the BASE data indices
ind_Ameriflux_Fluxnet = find(~isnan(t.AmeriFluxFLUXNETDataStart));
n_Ameriflux_Fluxnet = length(ind_Ameriflux_Fluxnet);

ind_Ameriflux_noFluxnet = find(isnan(t.AmeriFluxFLUXNETDataStart(ind_Ameriflux)));
n_Ameriflux_noFluxnet = length(ind_Ameriflux_noFluxnet);
% check
% n_Fluxnet_noFluxnet = n_Ameriflux_Fluxnet + n_Ameriflux_noFluxnet;



%% PLOTTING BEGINS HERE

%% . Canadian map of flux sites

% FIRST, run sections 1 and 2

addpath(genpath('/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux'));

saveplot = 0;   % 0 = no, 1 = yes

% load map outline
figure(1)
clf;
set(gcf,'color','white');


min_lat = 40;
max_lat = 85;
min_lon = -142;
% min_lon = -165;
max_lon = -54;
% max_lon = -45;

% load shapefile
% M = m_shaperead('geo_political_region_2');
M_nolakes = m_shaperead('ne_10m_admin_0_countries');
P = m_shaperead('ne_10m_admin_0_countries_lakes');

% this one for composite maps!
m_proj('lambert','longitudes',[min_lon max_lon], ...
           'latitudes',[min_lat max_lat]);%'direction','vertical','aspect',1);
% m_proj('miller','longitudes',[min_lon max_lon], ...
           % 'latitudes',[min_lat max_lat]);%'direction','vertical','aspect',1);

% define some colours
% RGB = [173 216 230]/256;
RGB_lightblue = [235 256 256]/256;
RGB_lightgrey = [240 240 240]/256;
RGB_grey = [211 211 211]/256;
RGB_darkgrey = [71 71 71]/256;

blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;
cl_colors = {blue, red, black, ...
             green, brown, purple};
cl_str_clr_names = ["blue", "red", "black", "green", "brown", "purple"];

RGB = orderedcolors("gem");
H = rgb2hex(RGB);

% plot coastline
m_coast('patch',RGB_grey,'lineWidth',0.5,'edgecolor',RGB_darkgrey);
% m_coast('patch',RGB_lightgrey,'lineWidth',0.5,'edgecolor',RGB_grey);
hold on
m_grid;

%------------------------------------------------------------------------
% plotting for CFI proposal
%------------------------------------------------------------------------

% find intersection for Ameriflux sites + active

% registered on Ameriflux + Fluxnet + INactive
common_inds = intersect(ind_Ameriflux_Fluxnet,ind_inactiveSites);
markerSize_other = 8;
for k = 1:length(common_inds)
   h6 = m_line(t.Longitude_degrees_(common_inds(k)),t.Latitude_degrees_(common_inds(k)),'marker','o','markersize',markerSize_other,'markeredgecolor','b','markerfacecolor','white','color',[0 0 0]);
end
% registered on Ameriflux + Fluxnet + active
common_inds = intersect(ind_Ameriflux_Fluxnet,ind_activeSites);
markerSize_other = 8;
for k = 1:length(common_inds)
   h5 = m_line(t.Longitude_degrees_(common_inds(k)),t.Latitude_degrees_(common_inds(k)),'marker','o','markersize',markerSize_other,'markeredgecolor','k','markerfacecolor','blue','color',[0 0 0]);
end

% registered on Ameriflux + NO Fluxnet + INactive
common_inds = intersect(ind_Ameriflux_noFluxnet,ind_inactiveSites);
markerSize_other = 8;
for k = 1:length(common_inds)
   h4 = m_line(t.Longitude_degrees_(common_inds(k)),t.Latitude_degrees_(common_inds(k)),'marker','square','markersize',markerSize_other,'markeredgecolor','g','markerfacecolor','white','color',[0 0 0]);
end
% registered on Ameriflux + NO Fluxnet + active
common_inds = intersect(ind_Ameriflux_noFluxnet,ind_activeSites);
markerSize_other = 8;
for k = 1:length(common_inds)
   h3 = m_line(t.Longitude_degrees_(common_inds(k)),t.Latitude_degrees_(common_inds(k)),'marker','square','markersize',markerSize_other,'markeredgecolor','k','markerfacecolor','green','color',[0 0 0]);
end

% NOT registered on Ameriflux + INactive
common_inds = intersect(ind_nonAmeriflux,ind_inactiveSites);
markerSize_other = 8;
for k = 1:length(common_inds)
   h4 = m_line(t.Longitude_degrees_(common_inds(k)),t.Latitude_degrees_(common_inds(k)),'marker','^','markersize',markerSize_other,'markeredgecolor','#FF8800','markerfacecolor','white','color',[0 0 0]);  % orange
end

% NOT registered on Ameriflux + active
common_inds = intersect(ind_nonAmeriflux,ind_activeSites);
markerSize_other = 8;
for k = 1:length(common_inds)
   h3 = m_line(t.Longitude_degrees_(common_inds(k)),t.Latitude_degrees_(common_inds(k)),'marker','^','markersize',markerSize_other,'markeredgecolor','k','markerfacecolor','#FF8800','color',[0 0 0]);  % orange
end

% remaining sites wth NO Ameriflux data
for k = 1:length(ind_noAmerifluxData_otherSites)
    h3_nodata = m_line(t.Longitude_degrees_(ind_otherSites(ind_noAmerifluxData_otherSites(k))), ...
                t.Latitude_degrees_(ind_otherSites(ind_noAmerifluxData_otherSites(k))), ...
                'marker','o','markersize',markerSize_other,'markeredgecolor','k','markerfacecolor','#DADADA','color',[0 0 0]);
end
% remaining sites WITH Ameriflux DATA
for k = 1:length(ind_AmerifluxData_otherSites)
    h3_data = m_line(t.Longitude_degrees_(ind_otherSites(ind_AmerifluxData_otherSites(k))), ...
                    t.Latitude_degrees_(ind_otherSites(ind_AmerifluxData_otherSites(k))), ...
                    'marker','o','markersize',markerSize_other,'markeredgecolor','k','markerfacecolor','#6a6a6a','color',[0 0 0]);
end


markerSize = 9;
% Ameriflux sites
for k = 1:ind-1
    h1 = m_line(t.Longitude_degrees_(k),t.Latitude_degrees_(k),'marker','^','markersize',markerSize,'markeredgecolor','k','markerfacecolor','#ff1493','color',[0 0 0]);
end

% non-Ameriflux sites
for k = ind+1:a
    h2 = m_line(t.Longitude_degrees_(k),t.Latitude_degrees_(k),'marker','o','markersize',markerSize,'markeredgecolor','k','markerfacecolor','#01c9fc','color',[0 0 0]);
end

% n_Ameriflux = length(find(~isnan(t.Latitude_degrees_(1:ind-1))));
% n_nonAmeriflux = length(find(~isnan(t.Latitude_degrees_(ind+1:a))));
% n_totalSites = n_Ameriflux + n_nonAmeriflux;

set(h1,'linestyle','none')
set(h2,'linestyle','none')
legend([h1 h2],{['Registered on Ameriflux (' num2str(n_Ameriflux) ')'],['Other Sites (' num2str(n_nonAmeriflux) ')']},'FontSize',20)

m_text(-40,65,[num2str(n_totalSites) ' sites with known coordinates'],'FontSize',16)








%% X. plot bar chart with site number information (Ameriflux/registered/BASE data)

[a,b] = size(t);

% total number of sites
numSites = length(find(~cellfun(@isempty, t.Name)));

% Ameriflux and non-Ameriflux sites
ind = find(strcmp(t.AmerifluxSiteID,'Sites likely not registered on Ameriflux'));   % identify index dividing Ameriflux and non-Ameriflux sites in spreadsheet
n_Ameriflux = length(find(~cellfun(@isempty,t.Name(1:ind-1))));
n_nonAmeriflux = length(find(~cellfun(@isempty,t.Name(ind+1:a))));

% sites with data
ind_notNan = ~isnan(t.NumberOfYearsOfAmeriFluxBASEData(1:ind-1));
ind_nonZero =  t.NumberOfYearsOfAmeriFluxBASEData(1:ind-1) > 0;
ind_withData = ind_notNan & ind_nonZero;
n_Ameriflux_withData = sum(ind_withData);
n_Ameriflux_noData = n_Ameriflux - n_Ameriflux_withData;

Categories = {'Total sites','Ameriflux (registered)','Ameriflux (data)','Ameriflux (no data)','Not registered'};
vals = [numSites; n_Ameriflux; n_Ameriflux_withData; n_Ameriflux_noData; n_nonAmeriflux];
% vals = [135;83;68;15;52];
percentVals = round((vals./numSites)*100);

figure(1);
clf;
set(gcf,'color','white');
% b = bar(vals',0.8);
b = bar(percentVals',0.8);
b.FaceColor = 'flat';
b.EdgeColor = 'none';

colors = [[1 31 75]./255;
            [3 57 108]./255;
            [0 91 150]./255;
            [100 151 177]./255;
            [179 205 224]./255];
b.CData = colors;

xtips = b.XEndPoints;
ytips = b.YEndPoints;
% labels = string(percentVals);
offset = 0.3;
for k = 1:length(vals)
    label = string(vals(k));
    label = strcat('(',label,')');
    % label = strcat(label,'%');
    text(xtips(k), ytips(k) + offset, label, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16);
end

set(gca,'YGrid','on')
set(gca,'XTickLabel',Categories,'FontSize',24,'XTickLabelRotation',36)

% ylim([0 145])
ax = gca;
ax.XAxis.TickLength = [0 0]; 
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 14;
ylabel('Frequency (%)','FontSize',24)
% ylabel('Frequency (sites)','FontSize',24)

