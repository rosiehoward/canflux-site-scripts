% Static map for CanFlux meeting
%
% Rosie Howard
% 2 October 2025

% Section list:
%   0. Extract data from relevant spreadsheet
%   1. Canadian map of flux sites (Ameriflux/non-Ameriflux)
%   2. Map of flux sites by ecosystem
%   3. Map of flux sites by flux measurement (CH4/CO2)
%   4. Bar chart of which fluxes (CO2/CH4) sample which ecosystems
%   5. Extract Canadian site information from Baldocchi 2025
%   6. Plot bar chart with site number information ('Total sites','Ameriflux (registered)','Ameriflux (data)','Ameriflux (no data)','Not registered')
%   7. Plot pie chart (donut) of ecosystems
%   8. Plot number of years of BASE data
%   9. Plot years of data (coupled with 8)
%   10. Plot/explore years of data per ecosystem
%   11. Economical analysis (very rough!)


%% 0. Extract data from TSV file
clear;
fpath = '/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux';
% fname = 'AmeriFlux-sites-CanadianORIGDOWNLOAD_19August2025.tsv';    % Ameriflux data
% fname_3Oct2025 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet3Oct2025_withLatLons.tsv';   % only sites with LatLons
% fname_7Oct2025 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet7Oct2025_edited.tsv';     % Ameriflux + self-reported
% fname_21Oct2025 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet21Oct2025.tsv';     % Ameriflux + self-reported
% fname_28Nov2025 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet28Nov2025.tsv';     % Ameriflux + self-reported
% fname_6Feb2026 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet6Feb2026.tsv';     % Ameriflux + self-reported
% fname_10Feb2026 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet10Feb2026_NotIncludingProposedSites.tsv';   % Ameriflux + self-reported (not including proposed sites) -> presented on Can-Peat 2026 poster
fname_20Feb2026 = 'Canadian-Flux-Sites_DownloadedFromGoogleSheet20Feb2026_NotIncludingProposedSites.tsv';   % Ameriflux + self-reported (not including proposed sites) -> for CFI proposal

fname_ORIG_17Sept2025 = 'Canadian-Flux-Sites_17Sept2025_ORIG_forComparison.tsv';


w = warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
% p = readtable(fullfile(fpath,fname), "FileType","text",'Delimiter', '\t');
t = readtable(fullfile(fpath,fname_20Feb2026), "FileType","text",'Delimiter', '\t');
v = readtable(fullfile(fpath,fname_ORIG_17Sept2025), "FileType","text",'Delimiter', '\t');

% fname_Baldocchi = 'CanadianSites_Baldocchi2025.csv';
% t = readtable(fullfile(fpath,fname_Baldocchi),'Delimiter', ',');

addpath(genpath('/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux'));

%% Identify CFI collaborators by comparing current site spreadsheet with original 
% (Canadian-Flux-Sites_17Sept2025_ORIG_forComparison), using lat/lon

% using latitudes
[out1, pos1] = ismember(t.Latitude_degrees_, v.Latitude_degrees_); 
% where pos == 0 --> new information = collaborator
isNew1 = (pos1 == 0);
collabs1 = t.PrincipalInvestigator__emailAddress_(isNew1);
collabs1 = unique(collabs1);

% using "Active site" column - i.e. has information provided (whether
% active = 1 or not = 0)
hasActive = ~isnan(t.ActiveSite);
collabs2 = t.PrincipalInvestigator__emailAddress_(hasActive);
collabs2 = unique(collabs2);

% find union of two collab lists
collabs_all = unique([collabs1(:); collabs2(:)], 'stable');

% write to file
writecell(collabs_all,fullfile(fpath,'CFI_siteOwnerCategories_contactList_fromSitesSpreadsheetComparison'),'QuoteStrings', 'none');

% using longitudes (gives the same collabs, can use latitudes)
% [out2, pos2] = ismember(t.Longitude_degrees_, v.Longitude_degrees_); 
% % where pos == 0 --> new information = collaborator
% isNew2 = (pos2 == 0);
% collabs2 = t.PrincipalInvestigator__emailAddress_(isNew2);
% collabs2 = unique(collabs2);

% using site names (this gives less new collabs, do not use)
% [out3, pos3] = ismember(t.Name, v.Name); 
% % where pos == 0 --> new information = collaborator
% isNew3 = (pos3 == 0);
% collabs3 = t.PrincipalInvestigator__emailAddress_(isNew3);
% collabs3 = unique(collabs3);

% using AmerifluxSiteID (does not work)
% [out4, pos4] = ismember(t.AmerifluxSiteID, v.AmerifluxSiteID); 
% % where pos == 0 --> new information = collaborator
% isNew4 = (pos4 == 0);
% collabs4 = t.PrincipalInvestigator__emailAddress_(isNew4);
% collabs4 = unique(collabs4);

% using FluxesMeasured (no better than latitude)
% [out5, pos5] = ismember(t.FluxesMeasured, v.FluxesMeasured); 
% % where pos == 0 --> new information = collaborator
% isNew5 = (pos5 == 0);
% collabs5 = t.PrincipalInvestigator__emailAddress_(isNew5);
% collabs5 = unique(collabs5);

%% load CFI contact categories

fname_contacts = 'CFI_siteOwnerCategories_contactList.tsv';
z = readtable(fullfile(fpath,fname_contacts), "FileType","text",'Delimiter', '\t');

% list of site PIs
sitePIs_CFI = z.sitePIs_CFI;
sitePIs_CFI(cellfun(@isempty,sitePIs_CFI))=[];
sitePIs_CFIemails = z.sitePIs_CFIemails;
sitePIs_CFIemails(cellfun(@isempty,sitePIs_CFIemails))=[];

% list of site collaborators
siteCollabs_CFI = z.siteCollabs_CFI;
siteCollabs_CFIemails = z.siteCollabs_CFIemails;


% find indices with sites owned by CFI PIs, and indices with sites owned by
% collaborators (survey responses who are not CFI PIs)
ind_CFI_PIs = find(contains(t.PrincipalInvestigator__emailAddress_,sitePIs_CFIemails,"IgnoreCase",true));
ind_collab_PIs = find(contains(t.PrincipalInvestigator__emailAddress_,siteCollabs_CFIemails,"IgnoreCase",true));

common_values = intersect(ind_CFI_PIs, ind_collab_PIs);
% for now it's Sean Carey (CFI PI) and Graham Clark (collab PI) both listed
% on four sites - remove these from collab list (keep on CFI PI list)
for i = 1:length(common_values)
    mask = (ind_collab_PIs ~= common_values(i));
    ind_collab_PIs = ind_collab_PIs(mask);
end

% find remaining indices - other sites not owned by CFI PIs or collabs
isNotInList = ~ismember(1:length(t.PrincipalInvestigator__emailAddress_), ind_CFI_PIs) & ~ismember(1:length(t.PrincipalInvestigator__emailAddress_), ind_collab_PIs);
ind_otherSites = find(isNotInList); 
% need to remove empty (where spreadsheet dividers are, Ameriflux/non-Ameriflux)
% this is the first index
% [a,b] = size(t);
% ind_divide = find(strcmp(t.AmerifluxSiteID,'Sites likely not registered on Ameriflux'));   % identify index dividing Ameriflux and non-Ameriflux sites in spreadsheet
% ind_nocoords = [1,ind_divide];

ind_nocoords = find(isnan(t.Latitude_degrees_));

% remove site listings if lat/lon doesn't exist (for all categories)
common_values = intersect(ind_CFI_PIs,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_CFI_PIs ~= common_values(i));
    ind_CFI_PIs = ind_CFI_PIs(mask);
end

common_values = intersect(ind_collab_PIs,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_collab_PIs ~= common_values(i));
    ind_collab_PIs = ind_collab_PIs(mask);
end

common_values = intersect(ind_otherSites,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_otherSites ~= common_values(i));
    ind_otherSites = ind_otherSites(mask);
end

%% highlight sites with data in Ameriflux (plotting in 1a below)

% find indices among CFI_PIs that have data
AmerifluxData_CFI_PIs = t.NumberOfYearsOfAmeriFluxBASEData(ind_CFI_PIs);
ind_AmerifluxData_CFI_PIs = find(AmerifluxData_CFI_PIs > 0);
% find indices among CFI_PIs that do NOT have data
N = length(ind_CFI_PIs); 
i1 = true(1, N); 
i1(ind_AmerifluxData_CFI_PIs) = false; 
ind_noAmerifluxData_CFI_PIs = find(i1); 

% find indices among collabs that have data
AmerifluxData_collabs = t.NumberOfYearsOfAmeriFluxBASEData(ind_collab_PIs);
ind_AmerifluxData_collabs = find(AmerifluxData_collabs > 0);
% find indices among collabs that do NOT have data
M = length(ind_collab_PIs); 
j1 = true(1, M); 
j1(ind_AmerifluxData_collabs) = false; 
ind_noAmerifluxData_collabs = find(j1); 

% find indices among remaining "other" sites that have data
AmerifluxData_otherSites = t.NumberOfYearsOfAmeriFluxBASEData(ind_otherSites);
ind_AmerifluxData_otherSites = find(AmerifluxData_otherSites > 0);
% find indices among remaining "other" sites that do NOT have data
P = length(ind_otherSites); 
k1 = true(1, P); 
k1(ind_AmerifluxData_otherSites) = false; 
ind_noAmerifluxData_otherSites = find(k1); 

% all sites, sites with data, new sites, unknown
n_totalSites = length(ind_CFI_PIs)+length(ind_collab_PIs)+length(ind_otherSites);
n_data = length(ind_AmerifluxData_CFI_PIs)+length(ind_AmerifluxData_collabs)+length(ind_AmerifluxData_otherSites);  % total sites with data
n_nodata = length(ind_noAmerifluxData_CFI_PIs)+length(ind_noAmerifluxData_collabs)+length(ind_noAmerifluxData_otherSites);  % total sites with NO data

%% active vs. historical vs. unknown (plotting in 1b below)

% ind_activeSites = (t.ActiveSite == 1);
% ind_inactiveSites = (t.ActiveSite == 0);
% ind_unknown = isnan(t.ActiveSite);
% 
% active_lats = t.Latitude_degrees_(ind_activeSites);
% active_lons = t.Longitude_degrees_(ind_activeSites);
% 
% inactive_lats = t.Latitude_degrees_(ind_inactiveSites);
% inactive_lons = t.Longitude_degrees_(ind_inactiveSites);
% 
% unknown_lats = t.Latitude_degrees_(ind_unknown);
% unknown_lons = t.Longitude_degrees_(ind_unknown);

ind_activeSites = find(t.ActiveSite == 1);
ind_inactiveSites = find(t.ActiveSite == 0);
ind_unknown = find(isnan(t.ActiveSite));

ind_nocoords = find(isnan(t.Latitude_degrees_));

% remove site listings if lat/lon doesn't exist (for all categories)
common_values = intersect(ind_activeSites,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_activeSites ~= common_values(i));
    ind_activeSites = ind_activeSites(mask);
end

% remove site listings if lat/lon doesn't exist (for all categories)
common_values = intersect(ind_inactiveSites,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_inactiveSites ~= common_values(i));
    ind_inactiveSites = ind_inactiveSites(mask);
end

% remove site listings if lat/lon doesn't exist (for all categories)
common_values = intersect(ind_unknown,ind_nocoords);
for i = 1:length(common_values)
    mask = (ind_unknown ~= common_values(i));
    ind_unknown = ind_unknown(mask);
end

%% PLOTTING BEGINS HERE

%% 1a. Canadian map of flux sites
% Does not need spreadsheet loading first (lat/lons saved separately)
% These will need updating as sites are added etc.

addpath(genpath('/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux'));

% load map outline
figure(1)
clf;
set(gcf,'color','white');


min_lat = 40;
max_lat = 85;
min_lon = -165;
max_lon = -45;

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
m_coast('patch',RGB_lightgrey,'lineWidth',0.5,'edgecolor',RGB_grey);
hold on
m_grid;

% for k=1:length(M.ncst)
%     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'linewidth',0.2,'color',RGB_grey);
% end
% 
% for k=1:length(M_nolakes.ncst)
%     m_line(M_nolakes.ncst{k}(:,1),M_nolakes.ncst{k}(:,2),'linewidth',0.2,'color','k');
% end

% *************************************
% % Ameriflux vs. non-Ameriflux sites
% *************************************
%
% % load site lat/lons
% load("BaldocchiLatLons.mat");
% % load("siteLatLons.mat");    % old - not sustainable, to many steps to update, read directly from downloaded table instead :)
% 
% [a,b] = size(t);
% ind = find(strcmp(t.AmerifluxSiteID,'Sites likely not registered on Ameriflux'));   % identify index dividing Ameriflux and non-Ameriflux sites in spreadsheet
% 
% markerSize = 9;
% % Ameriflux sites
% for k = 1:ind-1
%     h1 = m_line(t.Longitude_degrees_(k),t.Latitude_degrees_(k),'marker','^','markersize',markerSize,'markeredgecolor','k','markerfacecolor','#ff1493','color',[0 0 0]);
% end
% 
% % non-Ameriflux sites
% for k = ind+1:a
%     h2 = m_line(t.Longitude_degrees_(k),t.Latitude_degrees_(k),'marker','o','markersize',markerSize,'markeredgecolor','k','markerfacecolor','#01c9fc','color',[0 0 0]);
% end
% 
% n_Ameriflux = length(find(~isnan(t.Latitude_degrees_(1:ind-1))));
% n_nonAmeriflux = length(find(~isnan(t.Latitude_degrees_(ind+1:a))));
% 
% set(h1,'linestyle','none')
% set(h2,'linestyle','none')
% legend([h1 h2],{['Registered on Ameriflux (' num2str(n_Ameriflux) ')'],['Other Sites (' num2str(n_nonAmeriflux) ')']},'FontSize',20)

% plot "other" site first so they appear at back
markerSize_other = 8;
% % remaining sites (ALL)
% for k = 1:length(ind_otherSites)
%     h3 = m_line(t.Longitude_degrees_(ind_otherSites(k)),t.Latitude_degrees_(ind_otherSites(k)),'marker','o','markersize',markerSize_other,'markeredgecolor','k','markerfacecolor','#6a6a6a','color',[0 0 0]);
% end
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

% plot collabs next...
markerSize_collabs = 10;
% % collaborator sites (ALL)
% for k = 1:length(ind_collab_PIs)
%     h2 = m_line(t.Longitude_degrees_(ind_collab_PIs(k)),t.Latitude_degrees_(ind_collab_PIs(k)),'marker','v','markersize',markerSize_collabs,'markeredgecolor','k','markerfacecolor','#01c9fc','color',[0 0 0]);
% end
% collaborator sites with NO Ameriflux Data
for k = 1:length(ind_noAmerifluxData_collabs)
    h2_nodata = m_line(t.Longitude_degrees_(ind_collab_PIs(ind_noAmerifluxData_collabs(k))), ...
                t.Latitude_degrees_(ind_collab_PIs(ind_noAmerifluxData_collabs(k))), ...
                'marker','v','markersize',markerSize_collabs,'markeredgecolor','k','markerfacecolor','#CAE9F5','color',[0 0 0]);
end
% collaborator sites WITH AMERIFLUX DATA
for k = 1:length(ind_AmerifluxData_collabs)
    h2_data = m_line(t.Longitude_degrees_(ind_collab_PIs(ind_AmerifluxData_collabs(k))), ...
                t.Latitude_degrees_(ind_collab_PIs(ind_AmerifluxData_collabs(k))), ...
                'marker','v','markersize',markerSize_collabs,'markeredgecolor','k','markerfacecolor','#01c9fc','color',[0 0 0]);
end

% plot CFI PIs
markerSize_PIs = 12;
% % CFI PI sites (ALL)
% for k = 1:length(ind_CFI_PIs)
%     h1 = m_line(t.Longitude_degrees_(ind_CFI_PIs(k)),t.Latitude_degrees_(ind_CFI_PIs(k)),'marker','^','markersize',markerSize_PIs,'markeredgecolor','k','markerfacecolor','#ff1493','color',[0 0 0]);
% end
% CFI PI sites WITH NO AMERIFLUX DATA
for k = 1:length(ind_noAmerifluxData_CFI_PIs)
    h1_nodata = m_line(t.Longitude_degrees_(ind_CFI_PIs(ind_noAmerifluxData_CFI_PIs(k))), ...
                t.Latitude_degrees_(ind_CFI_PIs(ind_noAmerifluxData_CFI_PIs(k))), ...
                'marker','^','markersize',markerSize_PIs,'markeredgecolor','k','markerfacecolor','#FFE3F7','color',[0 0 0]);
end
% CFI PI sites WITH AMERIFLUX DATA
for k = 1:length(ind_AmerifluxData_CFI_PIs)
    h1_data = m_line(t.Longitude_degrees_(ind_CFI_PIs(ind_AmerifluxData_CFI_PIs(k))), ...
                t.Latitude_degrees_(ind_CFI_PIs(ind_AmerifluxData_CFI_PIs(k))), ...
                'marker','^','markersize',markerSize_PIs,'markeredgecolor','k','markerfacecolor','#ff1493','color',[0 0 0]);
end

%-----------------------------------------------
% legend and stats for plotting sites with data categories
%-----------------------------------------------
n_totalSites = length(ind_CFI_PIs)+length(ind_collab_PIs)+length(ind_otherSites);
n_data = length(ind_AmerifluxData_CFI_PIs)+length(ind_AmerifluxData_collabs)+length(ind_AmerifluxData_otherSites);  % total sites with data
n_nodata = length(ind_noAmerifluxData_CFI_PIs)+length(ind_noAmerifluxData_collabs)+length(ind_noAmerifluxData_otherSites);  % total sites with NO data

n_CFI_PIs = length(ind_CFI_PIs);    % total PI-owned    
% percent_CFI_PIs = round((n_CFI_PIs/n_totalSites)*100);  % percent PI-owned
n_CFI_PIs_data = length(ind_AmerifluxData_CFI_PIs);     % total PI-owned with data
percent_CFI_PIs_data = round((n_CFI_PIs_data/n_totalSites)*100);    % percentage of total sites that are PI-owned and have data
% percent_CFI_PIs_data_ofPIs = round((n_CFI_PIs_data/n_CFI_PIs_data)*100);    % percentage of PI-owned sites that have data
n_CFI_PIs_nodata = length(ind_noAmerifluxData_CFI_PIs);     % total PI-owned with NO data
percent_CFI_PIs_nodata = round((n_CFI_PIs_nodata/n_totalSites)*100);    % percentage of total sites that are PI-owned and have NO data
% percent_CFI_PIs_nodata_ofPIs = round((n_CFI_PIs_nodata/n_CFI_PIs_nodata)*100);    % percentage of PI-owned sites that have NO data

n_collab_PIs = length(ind_collab_PIs);  % total collab-owned
% percent_collab_PIs = round((n_collab_PIs/n_totalSites)*100);    % percent collab-owned
n_collab_PIs_data = length(ind_AmerifluxData_collabs);     % total PI-owned with data
percent_collab_PIs_data = round((n_collab_PIs_data/n_totalSites)*100);    % percentage of total sites that are PI-owned and have data
% percent_collab_PIs_data_ofcollabs = round((n_collab_PIs_data/n_collab_PIs_data)*100);    % percentage of PI-owned sites that have data
n_collab_PIs_nodata = length(ind_noAmerifluxData_collabs);     % total PI-owned with NO data
percent_collab_PIs_nodata = round((n_collab_PIs_nodata/n_totalSites)*100);    % percentage of total sites that are PI-owned and have NO data
% percent_collab_PIs_nodata_ofcollabs = round((n_collab_PIs_nodata/n_collab_PIs_nodata)*100);    % percentage of PI-owned sites that have NO data

n_otherSites = length(ind_otherSites);  % total other sites
% percent_otherSites = round((n_otherSites/n_totalSites)*100);   % percent other sites
n_otherSites_data = length(ind_AmerifluxData_otherSites);     % total PI-owned with data
percent_otherSites_data = round((n_otherSites_data/n_totalSites)*100);    % percentage of total sites that are PI-owned and have data
% percent_otherSites_data_ofother = round((n_otherSites_data/n_otherSites_data)*100);    % percentage of PI-owned sites that have data
n_otherSites_nodata = length(ind_noAmerifluxData_otherSites);     % total PI-owned with NO data
percent_otherSites_nodata = round((n_otherSites_nodata/n_totalSites)*100);    % percentage of total sites that are PI-owned and have NO data
% percent_otherSites_nodata_ofother = round((n_otherSites_nodata/n_otherSites_nodata)*100);    % percentage of PI-owned sites that have NO data

set(h1_data,'linestyle','none')
set(h1_nodata,'linestyle','none')
set(h2_data,'linestyle','none')
set(h2_nodata,'linestyle','none')
set(h3_data,'linestyle','none')
set(h3_nodata,'linestyle','none')
% set(d,'linestyle','none')
legend([h1_data h1_nodata h2_data h2_nodata h3_data h3_nodata],{ ...
    [num2str(n_CFI_PIs_data) ' CFI PI-owned with AmeriFlux data (' num2str(percent_CFI_PIs_data) '%)'], ...
    [num2str(n_CFI_PIs_nodata) ' CFI PI-owned no AmeriFlux data (' num2str(percent_CFI_PIs_nodata) '%)'], ...
    [num2str(n_collab_PIs_data) ' CFI collaborator-owned with AmeriFlux data (' num2str(percent_collab_PIs_data) '%)'], ...
    [num2str(n_collab_PIs_nodata) ' CFI collaborator-owned no AmeriFlux data (' num2str(percent_collab_PIs_nodata) '%)'], ...
    [num2str(n_otherSites_data) ' Other Sites with Ameriflux data (' num2str(percent_otherSites_data) '%)'], ...
    [num2str(n_otherSites_nodata) ' Other Sites no Ameriflux data (' num2str(percent_otherSites_nodata) '%)']},'FontSize',18);

m_text(-40,63,[num2str(n_totalSites) ' sites with known coordinates'],'FontSize',16)
m_text(-40,60,[num2str(n_data) ' sites with Ameriflux data'],'FontSize',16)


% legend([h3_nodata h3_data h2_nodata h2_data h1_nodata h1_data],{ ...
%     [num2str(n_otherSites_nodata) ' Other Sites no Ameriflux data (' num2str(percent_otherSites_nodata) '%)'], ...
%     [num2str(n_otherSites_data) ' Other Sites with Ameriflux data (' num2str(percent_otherSites_data) '%)'], ...
%     [num2str(n_collab_PIs_nodata) ' CFI collaborator-owned no AmeriFlux data (' num2str(percent_collab_PIs_nodata) '%)'], ...
%     [num2str(n_collab_PIs_data) ' CFI collaborator-owned with AmeriFlux data (' num2str(percent_collab_PIs_data) '%)'], ...
%     [num2str(n_CFI_PIs_nodata) ' CFI PI-owned no AmeriFlux data (' num2str(percent_CFI_PIs_nodata) '%)'], ...
%     [num2str(n_CFI_PIs_data) ' CFI PI-owned with AmeriFlux data (' num2str(percent_CFI_PIs_data) '%)']},'FontSize',18);


%---------------------------------------------------
% legend and stats for plotting all sites (no data categories)
%---------------------------------------------------
% n_totalSites = length(ind_CFI_PIs)+length(ind_collab_PIs)+length(ind_otherSites);
% n_CFI_PIs = length(ind_CFI_PIs);
% percent_CFI_PIs = round((n_CFI_PIs/n_totalSites)*100);
% n_collab_PIs = length(ind_collab_PIs);
% percent_collab_PIs = round((n_collab_PIs/n_totalSites)*100);
% n_otherSites = length(ind_otherSites);
% percent_otherSites = round((n_otherSites/n_totalSites)*100);
% 
% set(h1,'linestyle','none')
% set(h2,'linestyle','none')
% set(h3,'linestyle','none')
% % set(d,'linestyle','none')
% legend([h1 h2 h3],{[num2str(n_CFI_PIs) ' CFI PI-owned (' num2str(percent_CFI_PIs) '%)'], ...
%     [num2str(n_collab_PIs) ' CFI collaborator-owned (' num2str(percent_collab_PIs) '%)'], ...
%     [num2str(n_otherSites) ' Other Sites (' num2str(percent_otherSites) '%)']},'FontSize',20);
% 
% m_text(-32,65,[num2str(n_totalSites) ' sites with known coordinates'],'FontSize',16)



% Baldocchi data
% for k = 1:length(BaldocchiLatLons)
%     % m_line(siteLatLons(k,2),siteLatLons(k,1),'marker','^','markersize',8,'markeredgecolor','k','markerfacecolor','#ff1493');
%     h1 = m_line(BaldocchiLatLons(k,2),BaldocchiLatLons(k,1),'marker','^','markersize',8,'markeredgecolor','k','markerfacecolor','blue','color',[0 0 0]);
% end

% OLD
%Ameriflux/non-Ameriflux data
% for k = 1:length(siteLatLons_Ameriflux)
%     % m_line(siteLatLons(k,2),siteLatLons(k,1),'marker','^','markersize',8,'markeredgecolor','k','markerfacecolor','#ff1493');
%     h1 = m_line(siteLatLons_Ameriflux(k,2),siteLatLons_Ameriflux(k,1),'marker','^','markersize',8,'markeredgecolor','k','markerfacecolor','#ff1493','color',[0 0 0]);
% end
% 
% for k = 1:length(siteLatLons_nonAmeriflux)
%     h2 = m_line(siteLatLons_nonAmeriflux(k,2),siteLatLons_nonAmeriflux(k,1),'marker','o','markersize',8,'markeredgecolor','k','markerfacecolor','#01c9fc','color',[0 0 0]);
% end

%% 1b. Canadian map of flux sites - active vs. inactive

addpath(genpath('/Users/rosie/Documents/Micromet/CANFLUX_Database/Canadian_Flux_Sites/CanFlux'));

% load map outline
figure(1)
clf;
set(gcf,'color','white');

min_lat = 40;
max_lat = 85;
min_lon = -165;
max_lon = -45;

% load shapefile
% M = m_shaperead('geo_political_region_2');
M_nolakes = m_shaperead('ne_10m_admin_0_countries');
P = m_shaperead('ne_10m_admin_0_countries_lakes');
% m_proj('oblique mercator','lon',[min(lon_subset) max(lon_subset)],'lat',[min(lat_subset) max(lat_subset)],'direction','vertical','aspect',.5);

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
m_coast('patch',RGB_lightgrey,'lineWidth',0.5,'edgecolor',RGB_grey);
hold on
m_grid;

%****************************************
% Active vs. inactive vs. unknown sites *
%****************************************

markerSize = 8;
% Unknown sites (don't know if they are active or inactive)
for k = 1:length(ind_unknown)
    h1 = m_line(t.Longitude_degrees_(ind_unknown(k)),t.Latitude_degrees_(ind_unknown(k)),'marker','*','markersize',markerSize,'markeredgecolor','k','markerfacecolor','#6a6a6a','color',[0 0 0]);
end

markerSize = 10;
% inactive sites
for k = 1:length(ind_inactiveSites)
    h2 = m_line(t.Longitude_degrees_(ind_inactiveSites(k)),t.Latitude_degrees_(ind_inactiveSites(k)),'marker','square','markersize',markerSize,'markeredgecolor','k','markerfacecolor','#6F0000','color',[0 0 0]);
end

markerSize = 10;
% Active sites
for k = 1:length(ind_activeSites)
    h3 = m_line(t.Longitude_degrees_(ind_activeSites(k)),t.Latitude_degrees_(ind_activeSites(k)),'marker','o','markersize',markerSize,'markeredgecolor','k','markerfacecolor','#5EBB66','color',[0 0 0]);
end

%---------------------------------------------------
% legend and stats for plotting active/inactive/unknown
%---------------------------------------------------
n_totalSites = length(ind_activeSites)+length(ind_inactiveSites)+length(ind_unknown);
n_active = length(ind_activeSites);
percent_active = round((n_active/n_totalSites)*100);
n_inactive = length(ind_inactiveSites);
percent_inactive = round((n_inactive/n_totalSites)*100);
n_unknown = length(ind_unknown);
percent_unknown = round((n_unknown/n_totalSites)*100);

set(h1,'linestyle','none')
set(h2,'linestyle','none')
set(h3,'linestyle','none')
% set(d,'linestyle','none')
legend([h3 h2 h1],{[num2str(n_active) ' reported active (' num2str(percent_active) '%)'], ...
    [num2str(n_inactive) ' reported inactive (' num2str(percent_inactive) '%)'], ...
    [num2str(n_unknown) ' unknown (' num2str(percent_unknown) '%)'],},'FontSize',20);

m_text(-40,65,[num2str(n_totalSites) ' sites with known coordinates'],'FontSize',16)

%% 2. For mapping by ecosystem: alternative way to load from spreadsheet table

%**** Works with 28 Nov spreadsheet *********
% only for Ameriflux-defined sites currently...


% load map outline
figure(1)
clf;
set(gcf,'color','white');

min_lat = 40;
max_lat = 85;
min_lon = -165;
max_lon = -45;

% load shapefile
% M = m_shaperead('geo_political_region_2');
M_nolakes = m_shaperead('ne_10m_admin_0_countries');
P = m_shaperead('ne_10m_admin_0_countries_lakes');

m_proj('lambert','longitudes',[min_lon max_lon], ...
           'latitudes',[min_lat max_lat]);%'direction','vertical','aspect',1);

% define some colours
% RGB = [173 216 230]/256;
RGB_lightblue = [235 256 256]/256;
RGB_lightgrey = [240 240 240]/256;
RGB_grey = [211 211 211]/256;

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
m_coast('patch',RGB_lightgrey,'lineWidth',0.5,'edgecolor',RGB_grey);
hold on
m_grid;

% define data
SiteVegetation = t.VegetationDescription_IGBP_;    % Ameriflux + self-reported
classes = unique(SiteVegetation);
classes(cellfun(@isempty,classes)) = [];
percentVeg = NaN(length(classes),1);
classShort = cell(length(classes),1);
% extra = {'Water Bodies: Oceans, seas, lakes, reservoirs, and rivers. Can be either fresh or salt- water bodies.'};

for i = 1:length(classes)
    class = classes{i};
    numPerClass = length(find(strcmp(class,SiteVegetation)));
    percentVeg(i) = (numPerClass/length(SiteVegetation))*100;
    % get shortened class name
    tmp = split(classes{i},':');
    classShort{i} = strcat(tmp{1},[' (' num2str(numPerClass) ')']);
end

% percentVeg_rounded = round(percentVeg);

% T = table(percentVeg_rounded,classShort,mycolors);
% T_sort = sortrows(T,1,'descend');

for i = 1:length(classes)

    eval(['isEco' num2str(i) ' = find(strcmp(t.VegetationDescription_IGBP_,classes{' num2str(i) '}));']); %#ok<*EVLSEQVAR>
    eval(['Eco' num2str(i) '_LatLons = cat(2,t.Latitude_degrees_(isEco' num2str(i) '),t.Longitude_degrees_(isEco' num2str(i) '));']);

end
Ecos = {Eco1_LatLons,Eco2_LatLons,Eco3_LatLons,Eco4_LatLons,Eco5_LatLons,Eco6_LatLons,Eco7_LatLons,Eco8_LatLons,Eco9_LatLons};
mycolors = ["#8C1A66", "#D966B3", "#78AB30", "#C7E0A3", "#5B731F", "#B9E67D", "#954535", "#C19A6B",'#005B96']';
% myColors = ["#C7E0A3", "#C19A6B", "#D966B3", "#954535", "#5B731F", "#005B96", "#B9E67D", "#78AB30", "#8C1A66"]'; % sorted by number of sites
markers = {'<','pentagram','^','v','square','diamond','>','o','hexagram'};

% Ecosystem data
for j = 1:length(Ecos)
    for k = 1:length(Ecos{j})
        % h(j) = m_line(Ecos{j}(k,2),Ecos{j}(k,1),'marker',markers{j},'markersize',10,'markeredgecolor','k','markerfacecolor',mycolors{j},'color',[0 0 0]);
        h(j) = m_line(Ecos{j}(k,2),Ecos{j}(k,1),'marker','o','markersize',10,'markeredgecolor','k','markerfacecolor',mycolors{j},'color',[0 0 0]);
    end
end

set(h,'linestyle','none')
legend(h,classShort,'FontSize',14,'location','southwest')

%% 3. For mapping by flux measurements: CO2, CH4

%**** Works with 21 October spreadsheet *********

% load map outline
figure(1)
clf;
set(gcf,'color','white');

min_lat = 40;
max_lat = 85;
min_lon = -165;
max_lon = -45;

% load shapefile
% M = m_shaperead('geo_political_region_2');
M_nolakes = m_shaperead('ne_10m_admin_0_countries');
P = m_shaperead('ne_10m_admin_0_countries_lakes');

m_proj('lambert','longitudes',[min_lon max_lon], ...
           'latitudes',[min_lat max_lat]);%'direction','vertical','aspect',1);

% define some colours
% RGB = [173 216 230]/256;
RGB_lightblue = [235 256 256]/256;
RGB_lightgrey = [240 240 240]/256;
RGB_grey = [211 211 211]/256;

% plot coastline
m_coast('patch',RGB_lightgrey,'lineWidth',0.5,'edgecolor',RGB_grey);
hold on
m_grid;

% define data
allFluxes = t.FluxesMeasured;

tmp = [];
for i = 1:length(allFluxes)
    tmp1 = strrep(allFluxes{i},' ','');
    tmp2 = split(tmp1,',');
    % tmp3 = split(tmp2,':');
    tmp = [tmp; tmp2];
end
Fluxes = unique(tmp);



CO2_str = 'CO2';
FC_str = 'FC';
CH4_str = 'CH4';
FCH4_str = 'FCH4';
H2O_str = 'H2O';
O3_str = 'O3';

ind_Flux_Any = find(~cellfun(@isempty,t.FluxesMeasured));   % sites with any flux measured

% CO2
ind_1 = contains(t.FluxesMeasured,CO2_str);   % indices of sites with flux named by CO2
ind_2 = contains(t.FluxesMeasured,FC_str);     % indices of sites with flux named by FC
ind_CO2 = ind_1 | ind_2;                % all indices with either 'CO2' or 'FC' denoting the flux
numFC = sum(ind_CO2);                      % number of sites with CO2 measured

% CH4
ind_3 = contains(t.FluxesMeasured,CH4_str);   % indices of sites with flux named by CH4
ind_4 = contains(t.FluxesMeasured,FCH4_str); % indices of sites with flux named by FCH4
ind_CH4 = ind_3 | ind_4;              % all indices with either 'CH4' or 'FCH4' denoting the flux
numFCH4 = sum(ind_CH4);                    % number of sites with CH4 measured

% H2O
ind_H2O = contains(t.FluxesMeasured,H2O_str);   % indices of sites with flux named by FH2O or H2O
numFH2O = sum(ind_H2O);                         % number of sites with H2O measured

%O3
% so far all sites measuring O3 (=1) also measures FC, FCH4, H, LE
ind_O3 = contains(t.FluxesMeasured,O3_str);   % indices of sites with flux named by FO3 or O3
numFO3 = sum(ind_O3);                         % number of sites with O3 measured

% find actual indices for each flux
ind_array_CO2 = find(ind_CO2);             % array of actual indices for sites with CO2
ind_array_CH4 = find(ind_CH4);             % array of actual indices for sites with CH4
ind_array_H2O = find(ind_H2O);                    % array of actual indices for sites with H2O
ind_array_O3 = find(ind_O3);                    % array of actual indices for sites with O3

% find sites with ONLY CO2 and NOT CH4
ind_CO2_only = []; 
for i = 1:length(ind_array_CO2)
    ind_val = ind_array_CO2(i);
    if sum(ismember(ind_array_CH4,ind_val)) == 0
        ind_CO2_only(i) = ind_val;
    else
        ind_CO2_only(i) = NaN;
    end
end
numCO2_only = sum(~isnan(ind_CO2_only));
ind_CO2_only = ind_CO2_only(~isnan(ind_CO2_only));
ind_CO2_only = ind_CO2_only';

% so far all sites that measure CH4 also measure CO2
ind_Flux_CO2_CH4 = ind_CO2 & ind_CH4;
ind_array_CO2_CH4 = find(ind_Flux_CO2_CH4);
numFluxCO2_CH4 = sum(ind_Flux_CO2_CH4);

ind_Flux_CO2_H2O = ind_CO2 & ind_H2O;
ind_array_CO2_H2O = find(ind_Flux_CO2_H2O);
numFluxCO2_H2O = sum(ind_Flux_CO2_H2O);

ind_Flux_CO2_CH4_H2O = ind_CO2 & ind_CH4 & ind_H2O;
ind_array_CO2_CH4_H2O = find(ind_Flux_CO2_CH4_H2O);
numFluxCO2_CH4_H2O = sum(ind_Flux_CO2_CH4_H2O);

colors = [[1 31 75]./255;
            [3 57 108]./255;
            [0 91 150]./255;
            [100 151 177]./255;
            [179 205 224]./255];

% sites with FC
for k = 1:length(ind_CO2_only)
    ind = ind_CO2_only(k);
    h1 = m_line(t.Longitude_degrees_(ind),t.Latitude_degrees_(ind),'marker','^','markersize',9,'markeredgecolor','k','markerfacecolor',colors(3,:),'color',[0 0 0]);
end

% sites with FC and FCH4
for k = 1:length(ind_array_CO2_CH4)
    ind = ind_array_CO2_CH4(k);
    h2 = m_line(t.Longitude_degrees_(ind),t.Latitude_degrees_(ind),'marker','o','markersize',9,'markeredgecolor','k','markerfacecolor','#23e200','color',[0 0 0]);
end

set(h1,'linestyle','none')
set(h2,'linestyle','none')
legend([h1 h2],{['FC (' num2str(numCO2_only) ')'],['FC + FCH4 (' num2str(numFluxCO2_CH4) ')']},'FontSize',18)

%% 4. Investigate which fluxes (CO2/CH4) sample which ecosystems
% must have run section 3 first

% define data
% SiteVegetation = t.VegetationDescription_IGBP_;    % Ameriflux + self-reported
% classes = unique(SiteVegetation);
% classes(cellfun(@isempty,classes)) = [];

Ecos_CO2 = t.VegetationAbbreviation_IGBP_(ind_array_CO2);
allEcos = unique(Ecos_CO2);
allEcos(1) = [];    % not ideal way of doing this, works for now
allEcos(8) = [];    % not ideal way of doing this, works for now
Ecos_CO2_CH4 = t.VegetationAbbreviation_IGBP_(ind_array_CO2_CH4); 

Eco_CO2_only = sum(ind_CO2_only);
% Eco_CO2
Eco_Fluxes = [];
for k = 1:length(allEcos)
    class = allEcos{k};
    eval(['ind_' class '_CO2 = sum(contains(Ecos_CO2,class));']);
    eval(['ind_' class '_CO2_CH4 = sum(contains(Ecos_CO2_CH4,class));']);    
    eval(['Fluxes_' class ' = [ind_' class '_CO2, ind_' class '_CO2_CH4];']);
    eval(['Eco_Fluxes = [Eco_Fluxes; Fluxes_' class '];']);
end
Eco_Fluxes(:,3) = Eco_Fluxes(:,1) - Eco_Fluxes(:,2);

% plot Eco-Flux data
figure(4)
clf;
set(gcf,'color','white');
bar([Eco_Fluxes(:,3),Eco_Fluxes(:,2)],'stacked','EdgeColor','none')
set(gca,'XTickLabel',allEcos,'FontSize',14);%,'XTickLabelRotation',45)
set(gca,'YTick',0:2:36)
ylim([0 36])
set(gca,'YGrid','on')
legend(['FC (' num2str(numCO2_only) ')'],['FC + FCH4 (' num2str(numFluxCO2_CH4) ')'],'location','northwest')
xlabel('Ecosystem (Ameriflux)','FontSize',18)
ylabel('No. of sites','FontSize',18)

% colors = [[1 31 75]./255;
%             [3 57 108]./255;
%             [0 91 150]./255;
%             [100 151 177]./255;
%             [179 205 224]./255];
% b.CData = colors;

%% 5. Extract Canadian site information from Baldocchi 2025

CA_Sites_Orig = t.FluxnetCode;
CA_Sites = unique(CA_Sites_Orig);

% CA_Lat_Orig = t.Lat;
% CA_Lat = unique(CA_Lat_Orig);
% 
% CA_Lon_Orig = t.Lon;
% CA_Lon = unique(CA_Lon_Orig);

CA_Lat = NaN(size(CA_Sites));
CA_Lon = NaN(size(CA_Sites));

for i = 1:length(CA_Sites)
    ind = find(strcmp(CA_Sites{i},t.FluxnetCode));
    CA_Lat(i) = t.Lat(ind(1));  % take first matching index and get latitude
    CA_Lon(i) = t.Lon(ind(1));  % take first matching index and get latitude
end

T_Sites_Baldocchi = table(CA_Sites,CA_Lat,CA_Lon);
BaldocchiLatLons = table2array(T_Sites_Baldocchi(:,2:3));

% After removing duplicates and other symbols (e.g.,'??'), there are a
% total of 30 Canadian sites total



%% 6. plot bar chart with site number information

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



%% 7. Plot pie chart of ecosystems

% BASEVegetation = p.VegetationDescription_IGBP_; % Ameriflux only
SiteVegetation = t.VegetationDescription_IGBP_;    % Ameriflux + self-reported
% SiteVegetation = t.Var4;    % Ameriflux + self-reported
% SiteVegetation = t.VegetationAbbreviation_IGBP_;    % Ameriflux + self-reported

classes = unique(SiteVegetation);
% for i = 1:length(classes)
    % if length(classes{i}) > 3
        % classes{i} = {};
    % end
% end
classes(cellfun(@isempty,classes)) = [];
percentVeg = NaN(length(classes),1);
classShort = cell(length(classes),1);

for i = 1:length(classes)
    class = classes{i};
    numPerClass = length(find(strcmp(class,SiteVegetation)));
    percentVeg(i) = (numPerClass/length(SiteVegetation))*100;
    % get shortened class name
    tmp = split(classes{i},':');
    classShort{i} = strcat(tmp{1},[' (' num2str(numPerClass) ')']);
end

percentVeg_rounded = round(percentVeg);
mycolors = ["#8C1A66", "#D966B3", "#78AB30", "#C7E0A3", "#5B731F", "#B9E67D", "#954535", "#C19A6B",'#BEBEBE','#005B96']';
% mycolors = ["#8C1A66", "#D966B3", "#78AB30", "#C7E0A3", "#5B731F", "#B9E67D", "#954535", "#C19A6B",'#005B96']'; % old, does not include URB (urban)

T = table(percentVeg_rounded,classShort,mycolors);
T_sort = sortrows(T,1,'descend');

figure(2);
clf;
set(gcf,'color','white');
d = donutchart(T_sort.percentVeg_rounded,T_sort.classShort);
% d = donutchart(percentVeg_rounded,classShort);
% mycolors = ["#8C1A66", "#D966B3", "#78AB30", "#C7E0A3", "#5B731F", "#B9E67D", "#954535", "#C19A6B"];
% mycolorsOLD = ["#8C1A66", "#D966B3", "#78AB30", "#C7E0A3", "#5B731F", "#B9E67D", "#00A3A3", "#00E6E6"];
colororder(T_sort.mycolors);
% colororder(mycolors);
d.FontSize = 22;
d.EdgeColor = "white";
d.LineWidth = 4;
d.LabelStyle = "name";
d.CenterLabel = ["Canadian" "Flux Site" "Ecosystems"];
d.CenterLabelFontSize = 32;
d.StartAngle = -45;


%% 8. Plot number of years of BASE data
BASEData_NumYears = t.NumberOfYearsOfAmeriFluxBASEData;
BASEData_NumYears(BASEData_NumYears == 0) = NaN;

Data_moreThanFiveYrs = find(BASEData_NumYears > 5); % Altaf requested this information

figure(3);
% figure(1);
clf;
set(gcf,'color','white');
% subplot(1,5,1:2);
colors = [[1 31 75]./255;
            [3 57 108]./255;
            [0 91 150]./255;
            [100 151 177]./255;
            [179 205 224]./255];

nbins = 27;
% x = randn(10000,1);
[counts,edges] = histcounts(BASEData_NumYears,nbins);
% center = 0.5*(edges(1:end-1)+edges(2:end));
% b = bar(center, counts, 0.9);
b = bar(counts, 0.9);
b.FaceColor = colors(3,:);
% b.FaceColor = 'flat';
b.EdgeColor = 'none';

% h = histogram(BASEData_NumYears,'EdgeColor','none');
% h.EdgeColor = "none";         
xlim([0 28])
set(gca,'YGrid','on')
ax = gca;
ax.FontSize = 12;
ax.XAxis.TickLength = [0 0]; 

set(gca,'XTick',1:27)
% set(gca,'XTickLabel',{'',1:27})
% set(gca,'XTickLabel',{'','5','10','15','20','25'})

% xlabel('Number of years of BASE data','FontSize',18)
% xlabel('Length of BASE dataset (years)','FontSize',18)
ylabel('Number of sites','FontSize',18)
% ylabel('Frequency','FontSize',18)
title('Length of dataset in years','FontSize',18)
% title('Sites with BASE data by length of dataset','FontSize',18)


%% 9. Plot years of data

BASEData_Years = t.YearsOfAmeriFluxBASEData;
BASEData_Years_Vals = [];

count = 1;
for i = 1:length(BASEData_Years)
    tmp =  split(BASEData_Years{i,1},",");
    l_tmp = length(tmp);
    for j = 1:l_tmp
        if isempty(tmp{j})
            BASEData_Years_Vals(count) = NaN;
            count = count + 1;
        else
            BASEData_Years_Vals(count) = str2num(tmp{j});
            count = count + 1;
        end
    end

end

figure(4)
clf;
set(gcf,'color','white');

% subplot(1,5,3:5)

colors = [[1 31 75]./255;
            [3 57 108]./255;
            [0 91 150]./255;
            [100 151 177]./255;
            [179 205 224]./255];

nbins = 32;
[counts,edges] = histcounts(BASEData_Years_Vals);
% [counts,edges] = histcounts(BASEData_Years_Vals,nbins);
centre = 0.5*(edges(1:end-1)+edges(2:end));
b = bar(centre, counts, 0.9); 
% b.FaceColor = 'flat';
b.FaceColor = colors(3,:);
b.EdgeColor = 'none';
% histogram(BASEData_Years_Vals,'EdgeColor','none');
% h = histogram(BASEData_Years_Vals);
% h.BarWidth = 0.9;

xlim([1993 2026])
set(gca,'YGrid','on')
ax = gca;
ax.XAxis.TickLength = [0 0]; 
ax.FontSize = 12;
% xlabel('Year','FontSize',18)
ylabel('Number of sites','FontSize',18)
% ylabel('Frequency','FontSize',18)
title('Years with flux data','FontSize',18)
% title('Sites with BASE data by year','FontSize',18)

sgtitle('Canadian Flux Data on Ameriflux','FontSize',22)
% sgtitle('Historical Canadian Flux Data on Ameriflux','FontSize',22)?


%% 10. explore years of data per ecosystem
% load Ameriflux only for now (fname)

startYear = p.AmeriFluxBASEDataStart;
endYear = p.AmeriFluxBASEDataEnd;
totalNumYears = (endYear - startYear)+1;

rangeYears = p.YearsOfAmeriFluxBASEData;    % list of years directly from Ameriflux
BASEData_NumYears = p.NumberOfYearsOfAmeriFluxBASEData;
BASEData_NumYears(BASEData_NumYears == 0) = NaN;

% does range equal number of years? (i.e., could be year(s) missing
% somewhere within total range)
actualYears = {};
hypotheticalNumYears = {};
for i = 1:length(rangeYears)
    hypothetical = startYear(i):1:endYear(i);
    hypotheticalNumYears = [hypotheticalNumYears; hypothetical];

    actual = str2num(rangeYears{i});
    if isempty(actual)
        actualYears = [actualYears; NaN];
    else
        actualYears = [actualYears; actual];
    end
end
% now compare length of array for each year
for i = 1:length(rangeYears)
    diff(i) = abs(length(actualYears{i}) - length(hypotheticalNumYears{i}));
end

% after testing, can use start and end years and assume there is continuous
% data between those:
for i = 1:length(startYear)
    startDate(i) = datetime(startYear(i),1,1,"Format","dd-MMM-uuuu");
end
startDate = startDate';
for i = 1:length(endYear)
    endDate(i) = datetime(endYear(i),12,31,"Format","dd-MMM-uuuu");
end
endDate = endDate';

dates = cat(2,startDate,endDate);
P = max(cellfun(@length, actualYears));
allDataDates = NaN(length(rangeYears),P);

% data for plotting
T = table(p.SiteID,actualYears,BASEData_NumYears,p.VegetationAbbreviation_IGBP_);
T_sort = sortrows(T,3,'descend');
T_sort = rmmissing(T_sort);
[a,b] = size(T_sort);

buffer = 0.5;   % half a year for centering on whole years in plot

for i = 1:a
    years = T_sort.Var2{i};
    if isscalar(years)
        T_sort.Var5{i} = [years-buffer,years+buffer];
    else
        years(1) = years(1) - buffer;
        years(end) = years(end) + buffer;
        T_sort.Var5{i} = years(1):1:years(end);
    end
end

% % plot BASE data availability per site
% figure(1);
% clf;
% set(gcf,'color','white');
% for i = 1:a
%     plot(T_sort.Var5{i},i*ones(size(T_sort.Var5{i})),'LineWidth',10)
%     ylim([0 length(T_sort.Var5)+1])
%     xlim([1992 2026])
%     hold on
%     set(gca,'XGrid','on')
% end
% 
% set(gca,'XTick',1992:2:2026)
% set(gca,'YTick',1:1:a)
% set(gca,'YTickLabel',T_sort.Var1)
% 
% ylabel('Site ID','FontSize',16)
% xlabel('Year','FontSize',16)
% title('Ameriflux BASE Data Availability','FontSize',18)

% plot BASE data availability by ecosystem per site
% 'Evergreen Needleleaf Forests (33)'	"#C7E0A3"
% 'Permanent Wetlands (30)'	"#C19A6B"
% 'Croplands (12)'	"#D966B3"
% 'Open Shrublands (12)'	"#954535"
% 'Grasslands (6)'	"#5B731F"
% 'Water Bodies (5)'	"#005B96"
% 'Mixed Forests (4)'	"#B9E67D"
% 'Deciduous Broadleaf Forests (3)'	"#78AB30"
% 'Barren Sparse Vegetation (2)'	"#8C1A66"

T_sort_Eco = sortrows(T_sort,4,'ascend');

% separate into different tables
Ecos = unique(T_sort_Eco.Var4);
% mycolors = ["#C7E0A3", "#C19A6B", "#D966B3", "#954535", "#5B731F", "#B9E67D", "#78AB30", "#8C1A66"]';

mycolors = ['#8C1A66', "#D966B3", "#78AB30", "#C7E0A3", "#5B731F", "#B9E67D", "#954535", "#C19A6B"];

count = 1;
for i = 1:length(Ecos)
    ind = strcmp(T_sort_Eco.Var4, Ecos{i});
    T_sort_Eco{ind,'colors'} = mycolors(count);
    eval(['T_' Ecos{i} ' = T_sort_Eco(ind,:);']);
    count = count + 1;
end

figure(2);
clf;
set(gcf,'color','white');

[a,b] = size(T_sort_Eco);
for j = 1:length(Ecos)
    Eco = Ecos{j};
    for i = 1:a
        if strcmp(T_sort_Eco.Var4{i},Eco) == 0
            % plot(T_sort_Eco.Var5{i},i*ones(size(T_sort_Eco.Var5{i})),'LineWidth',10,'color',mycolors{j});
            continue
        else
            h(j) = plot(T_sort_Eco.Var5{i},i*ones(size(T_sort_Eco.Var5{i})),'LineWidth',10,'color',T_sort_Eco.colors{i});
        end
            hold on
    end
end

ylim([0 length(T_sort_Eco.Var5)+1])
xlim([1992 2026])
set(gca,'XGrid','on')
set(gca,'XTick',1992:1:2026)
set(gca,'YTick',1:1:a)
set(gca,'YTickLabel',T_sort_Eco.Var1)

ylabel('Site ID','FontSize',16)
xlabel('Year','FontSize',16)
title('Ameriflux BASE Data Availability by Ecosystem','FontSize',18)

% set(h,'linestyle','none')
l = legend(h,Ecos,'FontSize',14,'location','northwest');
l.Direction = 'Reverse';


%% 11. Economical analysis (very rough!)

% EC flux site approximate costs
cost_both_CO2_CH4_upper = 450000;        % approximate upper bound cost for EC site measuring both CO2 and CH4
cost_both_CO2_CH4_lower = 250000;        % approximate lower bound cost for EC site measuring both CO2 and CH4

cost_only_CO2_upper = 370000;            % approximate upper bound cost for EC site measuring only CO2
cost_only_CO2_lower = 190000;            % approximate lower bound cost for EC site measuring only CO2

% define data
allFluxes = t.FluxesMeasured;

% this part to identify which fluxes are listed (only needed once between
% spreadsheet changes)
tmp = [];
for i = 1:length(allFluxes)
    tmp1 = strrep(allFluxes{i},' ','');
    if contains(tmp1,':')
        tmp2 = split(tmp1,' ');
        tmp3 = split(tmp2,':');
        tmp = [tmp;tmp3];
    else
        tmp2 = split(tmp1,',');
        tmp = [tmp; tmp2];
    end
end
Fluxes = unique(tmp);

CO2_str = 'CO2';
FC_str = 'FC';
CH4_str = 'CH4';
FCH4_str = 'FCH4';
% H2O_str = 'H2O';
% O3_str = 'O3';

ind_Flux_Any = find(~cellfun(@isempty,t.FluxesMeasured));   % sites with any flux measured

% CO2
ind_1 = contains(t.FluxesMeasured,CO2_str);   % indices of sites with flux named by CO2
ind_2 = contains(t.FluxesMeasured,FC_str);     % indices of sites with flux named by FC
ind_CO2 = ind_1 | ind_2;                % all indices with either 'CO2' or 'FC' denoting the flux
numFC = sum(ind_CO2);                      % number of sites with CO2 measured

% CH4
ind_3 = contains(t.FluxesMeasured,CH4_str);   % indices of sites with flux named by CH4
ind_4 = contains(t.FluxesMeasured,FCH4_str); % indices of sites with flux named by FCH4
ind_CH4 = ind_3 | ind_4;              % all indices with either 'CH4' or 'FCH4' denoting the flux
numFCH4 = sum(ind_CH4);                    % number of sites with CH4 measured

% find actual indices for each flux
ind_array_CO2 = find(ind_CO2);             % array of actual indices for sites with CO2
ind_array_CH4 = find(ind_CH4);             % array of actual indices for sites with CH4

% find sites with ONLY CO2 and NOT CH4
ind_CO2_only = []; 
for i = 1:length(ind_array_CO2)
    ind_val = ind_array_CO2(i);
    if sum(ismember(ind_array_CH4,ind_val)) == 0
        ind_CO2_only(i) = ind_val;
    else
        ind_CO2_only(i) = NaN;
    end
end
numCO2_only = sum(~isnan(ind_CO2_only));
ind_CO2_only = ind_CO2_only(~isnan(ind_CO2_only));
ind_CO2_only = ind_CO2_only';

% so far all sites that measure CH4 also measure CO2
ind_Flux_CO2_CH4 = ind_CO2 & ind_CH4;
ind_array_CO2_CH4 = find(ind_Flux_CO2_CH4);
numFluxCO2_CH4 = sum(ind_Flux_CO2_CH4);

numSites = length(allFluxes) - 3;
numSitesUnknown = numSites - numFluxCO2_CH4 - numCO2_only;

% estimate cost ranges
% total upper bound = (num. CO2 only * cost of CO2_upper) + (num. CO2_CH4 * cost of CO2_CH4_upper) + numCO2_CH4 
% total lower bound = numCO2_only*costCO2_lower + numCO2_CH4*cost_CO2_CH4_lower

% upper total cost estimate for both CO2 and CH4 sites, assuming ratio of
% CO2+CH4:CO2_only is 1:3 (based on existing ratio)

ratio_CO2_CH4 = 41/113;     % "likelihood" of unknown sites being CO2 and CH4
ratio_CO2_only = 72/113;    % "likelihood" of unknown sites being CO2 only

total_cost_CO2_CH4_upper = cost_both_CO2_CH4_upper * numFluxCO2_CH4;      % total cost of all CO2+CH4 sites, at upper price
total_cost_CO2_CH4_lower = cost_both_CO2_CH4_lower * numFluxCO2_CH4;      % total cost of all CO2+CH4 sites, at lower price

total_cost_CO2_only_upper = cost_only_CO2_upper * numCO2_only;            % total cost of all CO2 only sites, at upper price
total_cost_CO2_only_lower = cost_only_CO2_lower * numCO2_only;            % total cost of all CO2 only sites at upper price

total_cost_unknown_CO2_CH4_upper = cost_both_CO2_CH4_upper * (numSitesUnknown*ratio_CO2_CH4);   % estimated cost of unknown sites likely to have CH4+CO2, at upper price
total_cost_unknown_CO2_CH4_lower = cost_both_CO2_CH4_lower * (numSitesUnknown*ratio_CO2_CH4);   % estimated cost of unknown sites likely to have CH4+CO2, at lower price

total_cost_unknown_CO2_only_upper = cost_only_CO2_upper * (numSitesUnknown*ratio_CO2_only);     % estimated cost of unknown sites likely to have CO2 only, at upper price
total_cost_unknown_CO2_only_lower = cost_only_CO2_lower * (numSitesUnknown*ratio_CO2_only);     % estimated cost of unknown sites likely to have CO2 only, at lower price

% total upper and lower costs
total_cost_upper = total_cost_CO2_CH4_upper + total_cost_CO2_only_upper ...
                    + total_cost_unknown_CO2_CH4_upper + total_cost_unknown_CO2_only_upper;
total_cost_upper = total_cost_upper*(1e-6);     % convert to millions

total_cost_lower = total_cost_CO2_CH4_lower + total_cost_CO2_only_lower ...
                    + total_cost_unknown_CO2_CH4_lower + total_cost_unknown_CO2_only_lower;
total_cost_lower = total_cost_lower*(1e-6);     % convert to millions

%% 12. plot categorical bar chart with PIs/collabs/other sites and data vs. no data

[a,b] = size(t);

% total number of sites
numSites = length(find(~cellfun(@isempty, t.Name)));

% Ameriflux and non-Ameriflux sites
% ind = find(strcmp(t.AmerifluxSiteID,'Sites likely not registered on Ameriflux'));   % identify index dividing Ameriflux and non-Ameriflux sites in spreadsheet
% n_Ameriflux = length(find(~cellfun(@isempty,t.Name(1:ind-1))));
% n_nonAmeriflux = length(find(~cellfun(@isempty,t.Name(ind+1:a))));
% 
% % sites with data
% ind_notNan = ~isnan(t.NumberOfYearsOfAmeriFluxBASEData(1:ind-1));
% ind_nonZero =  t.NumberOfYearsOfAmeriFluxBASEData(1:ind-1) > 0;
% ind_withData = ind_notNan & ind_nonZero;
% n_Ameriflux_withData = sum(ind_withData);
% n_Ameriflux_noData = n_Ameriflux - n_Ameriflux_withData;

Categories = {'Total sites','PI-owned (data)','PI-owned (no data)', ... 
                'Collab-owned (data)','Collab-owned (no data)', ... 
                'Other (data)','Other (no data)'};
vals = [numSites; n_CFI_PIs_data; n_CFI_PIs_nodata; n_collab_PIs_data; n_collab_PIs_nodata; n_otherSites_data; n_otherSites_nodata];
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
            [255 20 147]./255;
            [255 227 247]./255;
            [1 201 252]./255;
            [202 233 245]./255;
            [106 106 106]./255;
            [218 218 218]./255];
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



