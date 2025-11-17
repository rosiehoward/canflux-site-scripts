% createERA5_MSLPcomposites.m
% 
% Rosie Howard
% 14 December

% load data
clear;

% macbook or iMac
comp = 'rosiehoward';
% comp = 'rhoward';

addpath(['/Users/' comp '/Documents/UBC/WFRT/GrantsEtc/Elk_Valley_BCGov/obs_data/']);
addpath('/Volumes/WFRT-Ext26/era5/');

% pick cases, (c1) ERA5 or ERA5-Land, (c2) extreme or average
% c1 = 'ERA5-Land';
c1 = 'ERA5';
% c2 = 'average';
c2 = '95p';    % 95th percentile

% save data/plots?
save_data = 1;  % 0 = not now thank you, 1 = yes please
save_plot = 1;  % 0 = not now thank you, 1 = yes please

filepath = ['/Users/' comp '/Documents/UBC/WFRT/GrantsEtc/Elk_Valley_BCGov/obs_data/scripts/'];   % Macbook
% filepath = '/Users/rhoward/Documents/UBC/WFRT/GrantsEtc/Elk_Valley_BCGov/obs_data/scripts/';
savepath = ['/Users/' comp '/Documents/UBC/WFRT/GrantsEtc/Elk_Valley_BCGov/obs_data/ERA5_climatology/ERA5/figures/'];
datapath_clim = ['/Users/' comp '/Documents/UBC/WFRT/GrantsEtc/Elk_Valley_BCGov/obs_data/ERA5/ERA5_climatology/'];
% datapath = '/Volumes/WFRT-Ext26/era5/';
% datapath = '/Users/rhoward/Documents/UBC/WFRT/GrantsEtc/Elk_Valley_BCGov/obs_data/temp_data/';

months_leap = [31,29,31,30,31,30,31,31,30,31,30,31];
months_noleap = [31,28,31,30,31,30,31,31,30,31,30,31];
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};


%% open dates file to loop over

fid1=fopen([datapath_clim 'ERA5_clim_months_full.txt']);
% fid=fopen('ERA5_clim_dates.txt');
tline = fgetl(fid1);

 % initialize monthly array (should eventually be 61x101x33)
 monthly_MSLP_1 = [];    % MSLP monthly array
 monthly_MSLP_2 = [];    % MSLP monthly array
 monthly_MSLP_3 = [];    % MSLP monthly array
 monthly_MSLP_4 = [];    % MSLP monthly array
 monthly_MSLP_5 = [];    % MSLP monthly array
 monthly_MSLP_6 = [];    % MSLP monthly array
 monthly_MSLP_7 = [];    % MSLP monthly array
 monthly_MSLP_8 = [];    % MSLP monthly array
 monthly_MSLP_9 = [];    % MSLP monthly array
 monthly_MSLP_10 = [];    % MSLP monthly array
 monthly_MSLP_11 = [];    % MSLP monthly array
 monthly_MSLP_12 = [];    % MSLP monthly array

% loop over months
while ischar(tline)    % for month
   disp(tline)

    % determine if year is leap year
    year = str2double(tline(1:4));
    month = str2double(tline(5:6));
    if month < 10
        month_str = ['0' num2str(month)];
    else
        month_str = num2str(month);
    end

    if (mod(year,4) == 0 && ((mod(year,100) ~= 0) || (mod(year,400) == 0)))
        disp(['Working on year ' num2str(year) ', which is a leap year']);
        num_days = months_leap;
    else
        disp(['Working on year ' num2str(year) ', which is not a leap year']);
        num_days = months_noleap;
    end

    datapath1 = [datapath_clim 'composite_maps/monthlies/'];

    P = readmatrix([datapath1 'era5-' tline '_MSLP_avg.txt']);  % load data file for current year/month
       
    if month == 1
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    elseif month == 2
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    elseif month == 3
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    elseif month == 4
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    elseif month == 5
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    elseif month == 6
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']); 
    elseif month == 7
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);   
    elseif month == 8
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    elseif month == 9 
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']); 
    elseif month == 10
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    elseif month == 11 
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    elseif month == 12
        eval(['monthly_MSLP_' num2str(month) ' = cat(3,monthly_MSLP_' num2str(month) ', P);']);
    end
    
    clear P
    
    tline = fgetl(fid1); % read next file date
end

fclose(fid1);

% find average for each month over 33 years

for i = 1:length(months_leap)

    eval(['mean_monthly_MSLP_' num2str(i) ' = mean(monthly_MSLP_' num2str(i) ',3);']);
    
    % save monthly composite data
    if (save_data == 1)
        eval(['month_filetext = [''MSLP_avg_' num2str(i) '''];']); % for filename, if saving data
        eval(['writematrix(mean_monthly_MSLP_' num2str(i) ',[datapath1 month_filetext ''.txt''],''Delimiter'','','');']);
    end

    
end

mean_monthly_MSLP = cat(3,mean_monthly_MSLP_1,mean_monthly_MSLP_2,mean_monthly_MSLP_3,mean_monthly_MSLP_4,mean_monthly_MSLP_5, ...
                        mean_monthly_MSLP_6,mean_monthly_MSLP_7,mean_monthly_MSLP_8,mean_monthly_MSLP_9,mean_monthly_MSLP_10, ...
                        mean_monthly_MSLP_11,mean_monthly_MSLP_12);

%% create seasonal composites (means)

% loop over monthly files (not quite evenly weighted averages but close
% enough for now)
for i = 1:12
    eval(['MSLP_avg_' num2str(i) ' = readmatrix(''MSLP_avg_' num2str(i) '.txt'');']);      
end

MSLP_DJF_avg = cat(3,MSLP_avg_12,MSLP_avg_1,MSLP_avg_2);
MSLP_MAM_avg = cat(3,MSLP_avg_3,MSLP_avg_4,MSLP_avg_5);
MSLP_JJA_avg = cat(3,MSLP_avg_6,MSLP_avg_7,MSLP_avg_8);
MSLP_SON_avg = cat(3,MSLP_avg_9,MSLP_avg_10,MSLP_avg_11);

MSLP_DJF_avg = mean(MSLP_DJF_avg,3);
MSLP_MAM_avg = mean(MSLP_MAM_avg,3);
MSLP_JJA_avg = mean(MSLP_JJA_avg,3);
MSLP_SON_avg = mean(MSLP_SON_avg,3);

mean_seasonal_MSLP = cat(3,MSLP_DJF_avg,MSLP_MAM_avg,MSLP_JJA_avg,MSLP_SON_avg);

clear MSLP_avg_*

%% plot monthly or seasonal composites

% pick monthlies or seasonal
% datapath_comp = ['/Users/' comp '/Documents/UBC/WFRT/GrantsEtc/Elk_Valley_BCGov/obs_data/ERA5_climatology/composite_maps/monthlies/'];
% load([datapath_comp 'mean_monthly_MSLP.mat']);

datapath_comp = ['/Users/' comp '/Documents/UBC/WFRT/GrantsEtc/Elk_Valley_BCGov/obs_data/ERA5_climatology/ERA5/composite_maps/seasonal/'];
% load([datapath_comp 'mean_seasonal_MSLP.mat']);

res = 0.25; % ERA5 (ERA5-Land is 0.1)

lon = (0:res:359.75)';
lat = (-90:res:90)';

% select subset
inds_lat = 541:601;
inds_lon = 901:1001;

lon_subset = lon(inds_lon);
lon_subset = lon_subset - 360;
lat_subset = lat(inds_lat);
min_lon = min(lon_subset);
max_lon = max(lon_subset);
min_lat = min(lat_subset);
max_lat = max(lat_subset);
[LN,LT]=meshgrid(lon_subset,lat_subset);

% load shapefile
M = m_shaperead('geo_political_region_2');


% m_proj('oblique mercator','lon',[min(lon_subset) max(lon_subset)],'lat',[min(lat_subset) max(lat_subset)],'direction','vertical','aspect',.5);

% this one for composite maps!
m_proj('miller','longitudes',[min_lon max_lon], ...
           'latitudes',[min_lat max_lat]);%'direction','vertical','aspect',1);

% m_coast('linewidth',2,'color','k');

lim_low = min(min(min(mean_seasonal_MSLP)));
lim_high = max(max(max(mean_seasonal_MSLP)));
% lim_low = min(min(min(mean_monthly_MSLP)));
% lim_high = max(max(max(mean_monthly_MSLP)));
inc = 0.1;
levels = lim_low:inc:lim_high;
clim = [lim_low lim_high];

if (save_plot == 1)
    
close;
figure('units','centimeters','outerposition',[0 0 90 90]);      % seasonal
% figure('units','centimeters','outerposition',[0 0 90 45]);      % monthlies
set(gcf,'color','white');

% station coords
location_Pekisko = [-114.4156, 50.3683];
% location_MtTurnbull = [-114.8388889851049, 50.217778016813]; % old coords from Google Maps
% Met. station coords measured with GPS
location_MtTurnbull = [-114.8420598, 50.21780375];
location_CastleB = [-114.8268978, 50.15870688];
location_Chauncey = [-114.8035537, 50.13389935];
location_Ewin = [-114.7518536, 49.99548958];
location_Todhunter = [-114.779483, 50.106738];

int_method = 'linear';

cmap_u10 = cbrewer('div','RdBu',256,int_method);
cmap_v10 = cbrewer('div','PRGn',256,int_method);
cmap_d2m = cbrewer('seq','YlGnBu',256,int_method);
cmap_t2m = cbrewer('div','RdBu',256,int_method);
cmap_t2m = flipud(cmap_t2m);
cmap_ppn = cbrewer('seq','Blues',256,int_method);
cmap_z = cbrewer('seq','BuPu',256,int_method);
% cmap_z = cbrewer('div','BrBG',256,int_method);
% cmap_t2m = cbrewer('seq','YlGnBu','PCHIP');

% clim_t2m = [-20, 20];
%         ind_plot = 1:6:24;  % which time indices do you want to plot?
%         clim_t2m = [floor(min(min(min(t2m_C(:,:,ind_plot))))), ceil(max(max(max(t2m_C(:,:,ind_plot)))))];
%         clim_msl = [floor(min(min(min(msl_kPa(:,:,ind_plot))))), ceil(max(max(max(msl_kPa(:,:,ind_plot)))))];
% clim_sp = [floor(min(min(min(sp_kPa(:,:,ind_plot))))), ceil(max(max(max(sp_kPa(:,:,ind_plot)))))];
% clim_ppn = [0, ceil(max(max(max(tp_mm(:,:,ind_plot)))))];
%         clim_z = [0, ceil(max(max(z)))];    % truncate at "sea level" (m)
% clim_z = [floor(min(min(z))), ceil(max(max(z)))];

indexvalue_t2m = 0;             % centre of colormap - for background u-wind plot and w-wind

    % loop over seasons
    for i = 1:4
        
    % loop over months
%     for i = 1:length(months)
        
        % plot data
        var = squeeze(mean_seasonal_MSLP(:,:,i));
        var = flipud(var);
        %             h = image('XData',lon,'YData',lat,'CData',var2,'CDataMapping','scaled');     %#ok<NASGU>
        %             h = pcolor(lon(inds_lon), lat(inds_lat), var_subset);
        %     h = pcolor(lon, lat, z_t);
        %                 set(h,'edgeColor','none')
        
        subplot(2,2,i);     % seasonal
%         subplot(3,4,i);   % monthlies
        m_contourf(lon_subset,lat_subset,var,levels,'edgecolor','none');
        hold on;
%         grid on;
        daspect([1 1 1]);
        
        %             CustomMap1 = create_colormap(clim_z,indexvalue_t2m,cmap_z);
        CustomMap1 = create_colormap(clim,indexvalue_t2m,cmap_t2m); %,100,1,f_colorbar);
        colormap(CustomMap1);
        caxis(clim);
        h = colorbar;
        set(get(h,'label'),'string','MSLP (kPa)');
        h.FontSize = 20;
        m_grid;
%         set(gca,'XTick',225:5:250);
%         set(gca,'XTickLabel',-135:5:-110);
%         set(gca,'YTick',45:5:60);
%         set(gca,'YTickLabel',45:5:60);
        
        % draw coastline and boundaries
        for k=1:length(M.ncst)
            m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'linewidth',1,'color','k');
        end
        
         % seasonal
        if i == 1 || i == 3
            ylabel('Latitude','FontSize',22);
        end
        if i == 3 || i == 4
            xlabel('Longitude','FontSize',22);
        end
        
        % monthlies
%         if i == 1 || i == 5 || i == 9
%             ylabel('lat','FontSize',16);
%         end
%         if i == 9 || i == 10 || i == 11 || i == 12
%             xlabel('lon','FontSize',16);
%         end
        
        m_line(location_Pekisko(1),location_Pekisko(2),'marker','.','markersize',30,'color','k');
        m_line(location_Pekisko(1),location_Pekisko(2),'marker','.','markersize',10,'color','y');
        m_line(location_MtTurnbull(1),location_MtTurnbull(2),'marker','.','markersize',30,'color','k');
        m_line(location_MtTurnbull(1),location_MtTurnbull(2),'marker','.','markersize',10,'color','m');
        % m_line(location_CastleB(1),location_CastleB(2),'marker','.','markersize',12,'color','k');
        % m_line(location_CastleB(1),location_CastleB(2),'marker','.','markersize',4,'color','m');
        % m_line(location_Chauncey(1),location_Chauncey(2),'marker','.','markersize',12,'color','k');
        % m_line(location_Chauncey(1),location_Chauncey(2),'marker','.','markersize',4,'color','m');
        % m_line(location_Ewin(1),location_Ewin(2),'marker','.','markersize',12,'color','k');
        % m_line(location_Ewin(1),location_Ewin(2),'marker','.','markersize',4,'color','m');
        % m_line(location_Todhunter(1),location_Todhunter(2),'marker','.','markersize',12,'color','k');
        % m_line(location_Todhunter(1),location_Todhunter(2),'marker','.','markersize',4,'color','m');    


%         titletext =  months{i};   % monthlies

        if i == 1
            titletext = 'DJF';
        elseif i == 2
            titletext = 'MAM';
        elseif i == 3
            titletext = 'JJA';
        elseif i == 4
            titletext = 'SON';
        end

        title(titletext,'FontSize',24);
%         eval(['month_filetext = [''MSLP_avg_' num2str(i) '''];']); % for filename, if saving data
        
      
        % save plot
%         type = 'png';
%         im_res = 350;
%         str = ['print -d' type ' -r' num2str(im_res) ' ' datapath_comp month_filetext '.' type];
%         eval(str);
        
        
    end

%     sgtitletext = 'ERA5 MSLP monthly average (1990-2022)';    % monthlies
    % sgtitletext = 'ERA5 MSLP seasonal average (1990-2022)';     % seasonal
    % sgtitle(sgtitletext,'fontsize',2);
   
    % save plot
    filetext = 'MSLP_seasonal_composites';   % seasonal
%     filetext = 'MSLP_monthly_composites'; % monthlies
    type = 'png';
    im_res = 350;
    str = ['print -d' type ' -r' num2str(im_res) ' ' datapath_comp filetext '.' type];
    eval(str);

    
end




% EOF
