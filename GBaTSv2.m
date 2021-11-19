% GBaTSv2 Greenland likely basal thermal state version 2 analysis.
% 
% Joe MacGregor (NASA/GSFC)
% Last updated: 19 November 2021

clear

disp('Start with data loading...')

% switches or fixed values
decim                       = 5; % grid index decimation, indices but also effectively km
filt_type                   = 'triang'; % triang or exp, triang recommended by McCormack et al. (2019)
thick_filt_ref              = 10; % number of ice thicknesses over which to run filter
plotting                    = false;

% load various variables
load /Users/jamacgre/Documents/research/matlab/greenland/mat/grl_coast.mat num_coast x_coast y_coast % NSIDC Greenland coastline
load /Users/jamacgre/Documents/research/matlab/greenland/mat/core_int.mat lat_core lon_core % compiled Greenland ice-core locations
GBaTSv1                     = load('/Users/jamacgre/Documents/data/mat/bed.mat', 'frac_frozen_filled', 'frac_thawed_filled', 'frac_uncert_filled', 'mask_agree_searise_update', 'mask_likely_filled'); % GBaTSv1 datasets needed

% basal melt rate inferred from Greenland radiostratigraphy v1
RAD                         = load('/Users/jamacgre/Documents/research/matlab/greenland/mat/strain_all_9ka_krige.mat', 'melt_bed_grd', 'melt_bed_lo_grd', 'melt_bed_hi_grd'); % basal melt rate grids
load /Users/jamacgre/Documents/research/matlab/greenland/mat/mask_D1.mat mask_D1 x_D1_ord y_D1_ord % region of applicability of 1-D steady-state radiostratigraphic modeling (D<1)

% BedMachine v3, needed for SeaRISE to ISMIP6 comparison
BM3                         = struct;
BM3.x                       = 1e-3 .* double(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2017-09-20.nc', 'x'))'; % projected x, km
BM3.y                       = 1e-3 .* double(flipud(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2017-09-20.nc', 'y'))); % projected y, km
BM3.elev_bed                = rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2017-09-20.nc', 'bed')); % bed elevation, m
BM3.elev_bed_uncert         = rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2017-09-20.nc', 'errbed')); % bed elevation uncertainty, m
BM3.elev_surf               = double(rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2017-09-20.nc', 'surface'))); % surface elevation, m
BM3.mask_gris               = double(rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2017-09-20.nc', 'mask'))); % ice mask
BM3.thick                   = double(rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2017-09-20.nc', 'thickness'))); % ice thickness, m
BM3.thick(~BM3.thick)       = NaN; % NaN out zero ice thickness

% BedMachine v4
BM4                         = struct;
BM4.x                       = 1e-3 .* double(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2021-04-20.nc', 'x'))'; % projected x, km
BM4.y                       = 1e-3 .* double(flipud(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2021-04-20.nc', 'y'))); % projected y, km
BM4.elev_bed                = rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2021-04-20.nc', 'bed')); % bed elevation, m
BM4.elev_bed_uncert         = rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2021-04-20.nc', 'errbed')); % bed elevation uncertainty, m
BM4.elev_surf               = double(rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineG0reenland-2021-04-20.nc', 'surface'))); % surface elevation, m
BM4.mask_gris               = double(rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2021-04-20.nc', 'mask'))); % ice mask
BM4.thick                   = double(rot90(ncread('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2021-04-20.nc', 'thickness'))); % ice thickness, m
BM4.thick(~BM4.thick)       = NaN; % NaN out zero ice thickness

% simplified masks for whole Greenland maps
BM4.mask_combo              = BM4.mask_gris;
BM4.mask_gris(BM4.mask_gris ~= 2) ...
                            = 0;
BM4.mask_gris(BM4.mask_gris == 2) ...
                            = 1;
BM4.mask_gris               = logical(BM4.mask_gris); % now mask_gris is for ice sheet only
BM4.mask_combo(BM4.mask_combo == 3) ...
                            = 0;
BM4.mask_combo(BM4.mask_combo == 4) ...
                            = 1; % mask_combo includes 0/ocean 1/land 2/ice

% NaN out surface elevation on ocean
[BM3.elev_surf(~BM4.mask_combo), BM4.elev_surf(~BM4.mask_combo)] ...
                            = deal(NaN);

% MEaSURES multi-year surface speed mosaic, v1
SM1                         = load('/Users/jamacgre/Documents/research/matlab/greenland/mat/speed_greenland_measures_v1.mat', 'x_grd', 'y_grd', 'speed', 'speed_uncert', 'speed_x', 'speed_y', 'speed_x_uncert', 'speed_y_uncert');
mask_gris_SM1               = logical(interp2(BM4.x, BM4.y, BM4.mask_gris, SM1.x_grd, SM1.y_grd, 'nearest'));
[SM1.speed(~mask_gris_SM1), SM1.speed_uncert(~mask_gris_SM1), SM1.speed_x(~mask_gris_SM1), SM1.speed_y(~mask_gris_SM1), SM1.speed_x_uncert(~mask_gris_SM1), SM1.speed_y_uncert(~mask_gris_SM1)] ...
                            = deal(NaN); % remove NaNs based on BM4 mask

% Mouginot et al. (2019) large-scale drainage basins
M19                         = shaperead('/Users/jamacgre/Documents/data/GBaTSv2/Greenland_Basins_PS_v1.4.2_simple.shp'); 
num_M19                     = length(M19);
for ii = 1:num_M19
    [M19(ii).X, M19(ii).Y]  = deal((1e-3 .* M19(ii).X), (1e-3 .* M19(ii).Y));
end

% X/Y limits for standard map display, km
[x_min, x_max, y_min, y_max]= deal(-632, 846, -3344, -670);

% reference 1-km X/Y grid, km
[x_grd, y_grd]              = meshgrid(x_min:x_max, (y_min:y_max)');
[num_y_RAD, num_x_RAD]      = size(x_grd);

% parallels/meridians for Greenland maps in EPSG:3413
if license('checkout', 'map_toolbox')
    [~, ~, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland((1e3 .* [x_min x_max]), (1e3 .* [y_min y_max]), 5, 10);
else
    load /Users/jamacgre/Documents/research/matlab/greenland/mat/greenland_graticule x_paral y_paral x_merid y_merid
end

% letters for plots
letters                     = 'a':'z';

%% Load existing/new borehole and subglacial lake datasets

disp('Load old and new Colgan data archeology, Jordan et al. (2018), Bowling et al. (2019), Panton and Karlsson (2015) and Leysinger-Vieli et al. (2018) data...')

% old Colgan borehole summary first
S                           = shaperead('/Users/jamacgre/Documents/research/data/greenland/basal_evidence/borehole_freezing.shp');
[x_borehole_frozen, y_borehole_frozen] ...
                            = deal([S.X], [S.Y]);
temp_borehole_frozen        = [-11.8 -1.3 -6.0 -11.4 -100 -6.6 -100 -100 -100]; % Table 1 of MacGregor et al. (2016, JGR), used in ISMIP6 plot to colorize boreholes in Fig. 2, -100 if value unknown but frozen

S                           = shaperead('/Users/jamacgre/Documents/research/data/greenland/basal_evidence/boreholes_temperate.shp');
[x_borehole_thawed, y_borehole_thawed] ...
                            = deal([S.X], [S.Y]);
S                           = shaperead('/Users/jamacgre/Documents/research/data/greenland/basal_evidence/subglacial_lakes.shp');
[x_lake, y_lake]            = deal([S.X], [S.Y]);

% Liam's UPS projection, converted to EPSG:3413
load /Users/jamacgre/Documents/research/matlab/greenland/mat/gimp_proj.mat gimp_proj % EPSG:3413 projection
ups_n                       = gimp_proj;
ups_n.ProjParm(1:2)         = [90 0];
clear gimp_proj
[lat_borehole_frozen, lon_borehole_frozen] ...
                            = projinv(ups_n, x_borehole_frozen, y_borehole_frozen);
[lat_borehole_thawed, lon_borehole_thawed] ...
                            = projinv(ups_n, x_borehole_thawed, y_borehole_thawed);
[lat_borehole_thawed, lon_borehole_thawed] ...
                            = deal([lat_borehole_thawed lat_core(3)], [lon_borehole_thawed lon_core(3)]); % add EastGRIP
[lat_lake, lon_lake]        = projinv(ups_n, x_lake, y_lake);
[x_lake, y_lake]            = projfwd(projcrs(3413), lat_lake, lon_lake);
[x_lake, y_lake]            = deal((1e-3 .* x_lake), (1e-3 .* y_lake));

% switch NEEM from frozen to thawed
[lat_borehole_thawed, lon_borehole_thawed] ...
                            = deal([lat_borehole_thawed lat_borehole_frozen(2)], [lon_borehole_thawed lon_borehole_frozen(2)]);
[lat_borehole_frozen, lon_borehole_frozen, temp_borehole_frozen] ...
                            = deal(lat_borehole_frozen([1 3:end]), lon_borehole_frozen([1 3:end]), temp_borehole_frozen([1 3:end]));

% lat/lon to EPSG:3413
[x_borehole_frozen, y_borehole_frozen] ...
                            = projfwd(projcrs(3413), lat_borehole_frozen, lon_borehole_frozen);
[x_borehole_frozen, y_borehole_frozen] ...
                            = deal((1e-3 .* x_borehole_frozen), (1e-3 .* y_borehole_frozen));
[x_borehole_thawed, y_borehole_thawed] ...
                            = projfwd(projcrs(3413), lat_borehole_thawed, lon_borehole_thawed);
[x_borehole_thawed, y_borehole_thawed] ...
                            = deal((1e-3 .* x_borehole_thawed), (1e-3 .* y_borehole_thawed));
                        
% load Colgan et al. (2021) additional boreholes
C21                         = readtable('/Users/jamacgre/Documents/research/data/greenland/Additional_Boreholes_For_Joe.xlsx');
[x_borehole_frozen, y_borehole_frozen, temp_borehole_frozen] ...
                            = deal([x_borehole_frozen (1e-3 .* C21.EASTING([1:3 end]))'], [y_borehole_frozen (1e-3 .* C21.NORTHING([1:3 end]))'], [temp_borehole_frozen ([-2.3 -1 -1.3 -1.7] + (273.15 - pmp(C21.THICKNESS([1:3 end])')))]);
[x_borehole_thawed, y_borehole_thawed] ...
                            = deal([x_borehole_thawed (1e-3 .* C21.EASTING(4:(end - 1)))'], [y_borehole_thawed (1e-3 .* C21.NORTHING(4:(end - 1)))']);
[x_borehole_thawed(end), y_borehole_thawed(end)] ...
                            = projfwd(projcrs(3413), C21.LATITUDE(end - 1), C21.LONGITUDE(end - 1));
[x_borehole_thawed(end), y_borehole_thawed(end)] ...
                            = deal((1e-3 * x_borehole_thawed(end)), (1e-3 * y_borehole_thawed(end)));
[lat_borehole_thawed_new, lon_borehole_thawed_new] ...
                            = deal(70.75, -50.09);
[x_borehole_thawed_new, y_borehole_thawed_new] ...
                            = projfwd(projcrs(3413), lat_borehole_thawed_new, lon_borehole_thawed_new);
[x_borehole_thawed, y_borehole_thawed] ...
                            = deal([x_borehole_thawed (1e-3 * x_borehole_thawed_new)], [y_borehole_thawed (1e-3 * y_borehole_thawed_new)]); % add in Law et al. (2021) Store R20

% Bowling et al. (2019, Nature Communications) subglacial lakes
B19                         = load('/Users/jamacgre/Documents/research/data/greenland/Bowling_2019_lakes.mat', 'B19');
B19                         = B19.B19;

[x_lake, y_lake]            = deal([x_lake (1e-3 .* B19.x(5:6))'], [y_lake (1e-3 .* B19.y(5:6))']); % only include Bowling et al. (2019) subglacial lakes not within original Colgan dataset
B19                         = B19(7:end, :); % remove subglacial lakes from dataset to focus on basal water later

% Jordan et al. (2018, The Cryosphere) basal water estimates
J18                         = load('/Users/jamacgre/Documents/research/matlab/greenland/mat/Jordan_2018_basal_water.mat', 'x', 'y', 'basal_water');

% Leysinger-Vieli et al. (2018, Nature Communications) basal plumes
LV18                        = struct;
[LV18.x_plume_s, LV18.x_plume_l, LV18.y_plume_s, LV18.y_plume_l] ...
                            = deal(ncread('/Users/jamacgre/Documents/research/data/greenland/FreezeOnIndex_MappedPlume_Data.nc', 'plumesmall_x'), ...
                                   ncread('/Users/jamacgre/Documents/research/data/greenland/FreezeOnIndex_MappedPlume_Data.nc', 'plumelarge_x'), ...
                                   ncread('/Users/jamacgre/Documents/research/data/greenland/FreezeOnIndex_MappedPlume_Data.nc', 'plumesmall_y'), ...
                                   ncread('/Users/jamacgre/Documents/research/data/greenland/FreezeOnIndex_MappedPlume_Data.nc', 'plumelarge_y'));

% load Panton and Karlsson (2015) UDRs
PK15                        = load('/Users/jamacgre/Documents/research/matlab/greenland/mat/UDRs_PantonKarlsson2015.mat', 'latUDRall', 'lonUDRall', 'UDR');
[x_UDR, y_UDR]              = projfwd(projcrs(3413), PK15.latUDRall, PK15.lonUDRall); % lat/lon to EPSG:3413
[x_UDR, y_UDR]              = deal((1e-3 .* x_UDR), (1e-3 .* y_UDR)); % m to km

%% TRIM/DECIMATE DATA THAT SHOULDN'T BE FILTERED

disp('Trimming and decimating non-model data...')

% original grid sizes
[num_y_BM4, num_x_BM4]      = deal(length(BM4.y), length(BM4.x)); 
[num_y_SM1, num_x_SM1]      = size(SM1.x_grd);

% decimated grid
length_grd                  = 1e3 * decim; % grid size, m
ind_x                       = find(~mod(x_grd(1, :), decim), 1):decim:find(~mod(x_grd(1, :), decim), 1, 'last'); % column indices for decimated grid in reference 1-km grid
ind_y                       = find(~mod(y_grd(:, 1), decim), 1):decim:find(~mod(y_grd(:, 1), decim), 1, 'last'); % row indices for decimated grid in reference 1-km grid
[x_decim, y_decim]          = deal(x_grd(ind_y, ind_x), y_grd(ind_y, ind_x)); % decimated x/y grid, based on 1-km original

% decimated grid size
[num_decim_y, num_decim_x]  = size(x_decim);

% indices of decimated grid within original grids
ind_x_BM4                   = interp1(BM4.x, 1:num_x_BM4, x_decim(1, :), 'nearest');
ind_y_BM4                   = interp1(BM4.y, (1:num_y_BM4)', y_decim(:, 1), 'nearest')';
ind_x_SM1                   = interp1(SM1.x_grd(1, :), 1:num_x_SM1, x_decim(1, :), 'nearest');
ind_y_SM1                   = interp1(SM1.y_grd(:, 1), (1:num_y_SM1)', y_decim(:, 1), 'nearest');
ind_x_RAD                   = interp1(x_grd(1, :), 1:num_x_RAD, x_decim(1, :), 'nearest');
ind_y_RAD                   = interp1(y_grd(:, 1), (1:num_y_RAD)', y_decim(:, 1), 'nearest');

% just decimate masks that shouldn't be filtered
[mask_gris_decim, mask_combo_decim, mask_D1_decim] ...
                            = deal(BM4.mask_gris(ind_y_BM4, ind_x_BM4), BM4.mask_combo(ind_y_BM4, ind_x_BM4), mask_D1(ind_y_RAD, ind_x_RAD));

%% BUILD THICKNESS-DEPENDENT FILTERS

disp('Building smoothing filters...')

dxy_BM4                     = diff(BM4.x(1:2)); % 150 m grid
dxy_SM1                     = SM1.x_grd(1, 2) - SM1.x_grd(1, 1);  % 250 m grid
dxy_RAD                     = x_grd(1, 2) - x_grd(1, 1); % 1000 m grid

% build smoothing filter for different grid sizes
filt_decim_BM4              = cell(1, round((max(BM4.thick, [], 'all') * thick_filt_ref) / (1e3 * dxy_BM4)));
filt_decim_SM1              = cell(1, round((max(BM4.thick, [], 'all') * thick_filt_ref) / (1e3 * dxy_SM1)));
filt_decim_RAD              = cell(1, round((max(BM4.thick, [], 'all') * thick_filt_ref) / (1e3 * dxy_RAD)));
[num_filt_BM4, num_filt_SM1, num_filt_RAD] ...
                            = deal(length(filt_decim_BM4), length(filt_decim_SM1), length(filt_decim_RAD));

% loop through to build filters, based on smallest grid size (BM4)
for ii = 1:num_filt_BM4
    [tmp1, tmp2]            = meshgrid([-ii:0 -1:-1:-ii], ([-ii:0 -1:-1:-ii])'); % temporary index i/j grid
    switch filt_type
        case 'triang'
            filt_decim_BM4{ii} ...
                            = interp1([ii 0], [0 1], sqrt((tmp1 .^ 2) + (tmp2 .^ 2))); % follow McCormack et al. (2019, Polar Res.)
        case 'exp'
            filt_decim_BM4{ii} ...
                            = interp1([0 ii], [0 1], exp(-sqrt((tmp1 .^ 2) + (tmp2 .^ 2)) ./ ii)); % 2-D thickness-dependent exponential filter
            filt_decim_BM4{ii}(~logical(interp1([ii 0], [0 1], sqrt((tmp1 .^ 2) + (tmp2 .^ 2)), 'linear', 0))) ...
                            = NaN; % NaN out values outside radius of interest
    end
    filt_decim_BM4{ii}      = filt_decim_BM4{ii} ./ sum(filt_decim_BM4{ii}, 'all', 'omitnan'); % normalize so that sum = 1
    if (ii <= num_filt_SM1)
        filt_decim_SM1{ii}  = filt_decim_BM4{ii};
    end
    if (ii <= num_filt_RAD)
        filt_decim_RAD{ii}  = filt_decim_BM4{ii};
    end
end

%% FILTER DATA

disp('Filtering selected grids...')

% filter input non-model decimated grids (vars) into filtered outputs (vars_filt)
vars                        = {'BM4.elev_bed_uncert'  'BM4.elev_surf'  'BM4.thick'  'SM1.speed_x'  'SM1.speed_y'  'SM1.speed_x_uncert'  'SM1.speed_y_uncert' 'RAD.melt_bed_grd' 'RAD.melt_bed_lo_grd' 'RAD.melt_bed_hi_grd'};  % input names of grids to be filtered
vars_filt                   = {'elev_bed_uncert_filt' 'elev_surf_filt' 'thick_filt' 'speed_x_filt' 'speed_y_filt' 'speed_x_uncert_filt' 'speed_y_uncert_filt' 'melt_bed_filt' 'melt_bed_lo_filt' 'melt_bed_hi_filt'}; % filtered output names
num_var                     = length(vars); % number of variables to filter

% initialize to NaN
for ii = 1:num_var
    eval([vars_filt{ii} ' = NaN(num_decim_y, num_decim_x);'])
end

% closest matching set of ice thicknesses to thicknesses at decimated grid locations
filt_decim_ref_BM4          = interp1(((1:length(filt_decim_BM4)) .* ((1e3 * dxy_BM4) / thick_filt_ref)), 1:length(filt_decim_BM4), BM4.thick(ind_y_BM4, ind_x_BM4), 'nearest', 'extrap');
filt_decim_ref_SM1          = interp1(((1:length(filt_decim_SM1)) .* ((1e3 * dxy_SM1) / thick_filt_ref)), 1:length(filt_decim_SM1), BM4.thick(ind_y_BM4, ind_x_BM4), 'nearest', 'extrap');
filt_decim_ref_RAD          = interp1(((1:length(filt_decim_RAD)) .* ((1e3 * dxy_RAD) / thick_filt_ref)), 1:length(filt_decim_RAD), BM4.thick(ind_y_BM4, ind_x_BM4), 'nearest', 'extrap');
[filt_decim_ref_BM4(~mask_gris_decim), filt_decim_ref_SM1(~mask_gris_decim), filt_decim_ref_RAD(~mask_gris_decim)] ...
                            = deal(NaN);

for jj = 1:num_decim_x % loop through all columns
    
    if ~mod(jj, 25)
        disp([num2str(round(1e2 * (jj / num_decim_x))) '%']) % progress indicator
    end
    
    for ii = find(mask_gris_decim(:, jj))' % only bother with columns in current row that pass the mask test
        
        % current filter matrix and half-size
        filt_decim_curr_BM4 = filt_decim_BM4{filt_decim_ref_BM4(ii, jj)};
        num_filt_curr_BM4   = (size(filt_decim_curr_BM4, 1) - 1) / 2;
        
        % indices within original matrix
        ind_tmp_BM4         = -num_filt_curr_BM4:num_filt_curr_BM4; % reference index set
        ii_tmp_BM4          = repmat((ind_y_BM4(ii) + ind_tmp_BM4)', ((2 * num_filt_curr_BM4) + 1), 1); % row indices to filter
        jj_tmp_BM4          = ind_x_BM4(jj) + ind_tmp_BM4(ones(((2 * num_filt_curr_BM4) + 1), 1), :); % column indices to filter
        jj_tmp_BM4          = jj_tmp_BM4(:); % column-ify
        ind_good_BM4        = find((ii_tmp_BM4 > 0) & (jj_tmp_BM4 > 0) & (ii_tmp_BM4 <= num_y_BM4) & (jj_tmp_BM4 <= num_x_BM4)); % deal with edges
        ind_curr_BM4        = sub2ind([num_y_BM4 num_x_BM4], ii_tmp_BM4(ind_good_BM4), jj_tmp_BM4(ind_good_BM4)); % convert row/column indices to linear indices
        
        % repeat for speed grids
        filt_decim_curr_SM1 = filt_decim_SM1{filt_decim_ref_SM1(ii, jj)};
        num_filt_curr_SM1   = (size(filt_decim_curr_SM1, 1) - 1) / 2;
        
        ind_tmp_SM1         = -num_filt_curr_SM1:num_filt_curr_SM1;
        ii_tmp_SM1          = repmat((ind_y_SM1(ii) + ind_tmp_SM1)', ((2 * num_filt_curr_SM1) + 1), 1);
        jj_tmp_SM1          = ind_x_SM1(jj) + ind_tmp_SM1(ones(((2 * num_filt_curr_SM1) + 1), 1), :);
        jj_tmp_SM1          = jj_tmp_SM1(:);
        ind_good_SM1        = find((ii_tmp_SM1 > 0) & (jj_tmp_SM1 > 0) & (ii_tmp_SM1 <= num_y_SM1) & (jj_tmp_SM1 <= num_x_SM1));
        ind_curr_SM1        = sub2ind([num_y_SM1 num_x_SM1], ii_tmp_SM1(ind_good_SM1), jj_tmp_SM1(ind_good_SM1));
        
        % repeat for radar grids
        filt_decim_curr_RAD = filt_decim_RAD{filt_decim_ref_RAD(ii, jj)};
        num_filt_curr_RAD   = (size(filt_decim_curr_RAD, 1) - 1) / 2;
        
        ind_tmp_RAD         = -num_filt_curr_RAD:num_filt_curr_RAD;
        ii_tmp_RAD          = repmat((ind_y_RAD(ii) + ind_tmp_RAD)', ((2 * num_filt_curr_RAD) + 1), 1);
        jj_tmp_RAD          = ind_x_RAD(jj) + ind_tmp_RAD(ones(((2 * num_filt_curr_RAD) + 1), 1), :);
        jj_tmp_RAD          = jj_tmp_RAD(:);
        ind_good_RAD        = find((ii_tmp_RAD > 0) & (jj_tmp_RAD > 0) & (ii_tmp_RAD <= num_y_RAD) & (jj_tmp_RAD <= num_x_RAD));
        ind_curr_RAD        = sub2ind([num_y_RAD num_x_RAD], ii_tmp_RAD(ind_good_RAD), jj_tmp_RAD(ind_good_RAD));
        
        % filter each variable by multiplying thickness-appropriate normalized 2-D matrix with local values and summing result 
        for kk = 1:num_var
            if (eval(['sum(filt_decim_curr_' vars{kk}(1:3) '(~isnan(' vars{kk} '(ind_curr_' vars{kk}(1:3) '))), ''omitnan'');']) >= 0.5) % limits bias at edge of grid
                eval([vars_filt{kk} '(ii, jj) = sum((filt_decim_curr_' vars{kk}(1:3) '(ind_good_' vars{kk}(1:3) ') .* ' vars{kk} '(ind_curr_' vars{kk}(1:3) ')), ''omitnan'') / sum_filt_decim_curr;'])
            end
        end
    end
end

%% SPEED FILTERING VIA AZIMUTH AND DEFORMATION SPEED CALCULATION

disp('Filter surface speed with elevation azimuth and calculate deformation speed...')

% filtered modern surface speed and uncertainty
speed_filt                  = sqrt((speed_x_filt .^ 2) + (speed_y_filt .^ 2)); % m/yr
speed_uncert_filt           = sqrt((speed_x_uncert_filt .^ 2) + (speed_y_uncert_filt .^ 2));

% min/max speeds
speed_filt_min              = speed_filt - speed_uncert_filt;
speed_filt_min(speed_filt_min <= 0) ...
                            = 0;
speed_filt_max              = speed_filt + speed_uncert_filt;

% apply GrIS mask to filtered surface elevation
elev_surf_filt(~mask_gris_decim) ...
                            = NaN;

% surface-velocity and elevation-gradient flow azimuths
az_speed                    = atan2(speed_y_filt, speed_x_filt); % rad
[elev_grad_x, elev_grad_y]  = gradient(elev_surf_filt, length_grd, length_grd);
az_elev                     = atan2(-elev_grad_y, -elev_grad_x); % rad

% extract trigonometric elements
[az_sin_cat, az_cos_cat]    = deal(NaN(num_decim_y, num_decim_x, 2));
az_sin_cat(:, :, 1)         = sin(az_speed);
az_cos_cat(:, :, 1)         = cos(az_speed);
az_sin_cat(:, :, 2)         = sin(az_elev);
az_cos_cat(:, :, 2)         = cos(az_elev);

% weight filter exponentially using reference speed
speed_az_decay              = 100; % speed above which weight is unity for InSAR surface speeds, m/s (100 in GBTSv1)
speed_uncert_rel_decay      = 0.25; % (0.1 in GBTSv1, 0.25 here), also fractional cutoff for speeds
wt_az_elev                  = exp(-speed_filt ./ speed_az_decay);% + exp(-speed_uncert_rel_decay ./ (speed_uncert_filt ./ speed_filt)); % characteristic length of surface speed to unity ratio
wt_az_elev(wt_az_elev > 1)  = 1;
wt_az_speed                 = 1 - wt_az_elev;
wt_az_elev(isnan(az_speed)) = 1; % maximum weight (1) if the other is NaN
wt_az_speed(isnan(az_elev)) = 1;
az_mean                     = atan2(((az_sin_cat(:, :, 1) .* wt_az_speed) + (az_sin_cat(:, :, 2) .* wt_az_elev)), ((az_cos_cat(:, :, 1) .* wt_az_speed) + (az_cos_cat(:, :, 2) .* wt_az_elev))); % mean azimuth, radians
az_mean_cos                 = cos(az_mean); % rad
az_mean_sin                 = sin(az_mean); % rad

% % reproject speed in x/y directions (not currently used, but nice to have)
% [speed_x_filt2, speed_y_filt2] ...
%                             = deal((speed_filt .* az_mean_cos), (speed_filt .* az_mean_sin)); % reprojected x/y filtered speeds and uncertainties

% deformation speed calculations
density_ice                 = 900; % ice column density, kg m^-3
gravity_std                 = 9.81; % acceleration due to gravity at Earth's surface, m/(s^2)
slope_surf_flow             = abs((elev_grad_x .* az_mean_cos) + (elev_grad_y .* az_mean_sin)); % project surface slope onto weighted flow direction
driving_stress              = (density_ice * gravity_std) .* thick_filt .* slope_surf_flow; % standard driving stress for sloped slab of ice along flow, Pa
driving_stress_lo           = (density_ice * gravity_std) .* (thick_filt - elev_bed_uncert_filt) .* slope_surf_flow; % lo value
driving_stress_hi           = (density_ice * gravity_std) .* (thick_filt + elev_bed_uncert_filt) .* slope_surf_flow; % hi value
E                           = 2; % enhancement factor, dimensionless
E_min                       = 1;
E_max                       = 4;
n                           = 3; % default Glen's flow law rate exponent, dimensionless
speed_def                   = ((2 * E * ice_visc(273.15)) / (n + 1)) .* (driving_stress .^ n) .* thick_filt; % deformation speed, p. 310 in Cuffey and Paterson (2010)
% speed_def_uncert            = speed_def .* sqrt(((4 .* (elev_bed_uncert_filt ./ thick_filt)) .^ 2) + (E_uncert_rel .^ 2)); % propagation of error, no longer used because of non-monotonic range of E
speed_def_lo                = ((2 * E_min * ice_visc(273.15)) / (n + 1)) .* (driving_stress .^ n) .* (thick_filt - elev_bed_uncert_filt);
speed_def_lo(speed_def_lo < 0) ...
                            = 0;
speed_def_hi                = ((2 * E_max * ice_visc(273.15)) / (n + 1)) .* (driving_stress .^ n) .* (thick_filt + elev_bed_uncert_filt);
speed_def_n4                = ((2 * E * (3.3e-29 * (60 * 60 * 24 * 365.25))) / (4 + 1)) .* (driving_stress .^ 4) .* thick_filt; % Bons et al. (2018) deformation speed, m/yr using

% surface speed over deformation speed ratios
speed_ratio_std             = speed_filt ./ speed_def; % std=standard output
speed_ratio_min             = speed_filt_min ./ speed_def_hi; % min ratio is min numerator over max denominator
speed_ratio_max             = speed_filt_max ./ speed_def_lo; % max ratio is max numerator over min denominator

% indices where surface speed uncertainty above threshold
ind_speed_uncert_cutoff     = find((speed_uncert_filt ./ speed_filt) >= speed_uncert_rel_decay);

% trim ratios where surface speed uncertainty renders interpretation dubious
[speed_ratio_std_trim, speed_ratio_min_trim, speed_ratio_max_trim] ...
                            = deal(speed_ratio_std, speed_ratio_min, speed_ratio_max);
[speed_ratio_std_trim(ind_speed_uncert_cutoff), speed_ratio_min_trim(ind_speed_uncert_cutoff), speed_ratio_max_trim(ind_speed_uncert_cutoff)] ...
                            = deal(NaN);
[speed_ratio_std_trim(speed_ratio_std_trim > 20), speed_ratio_min_trim(speed_ratio_min_trim > 20), speed_ratio_max_trim(speed_ratio_max_trim > 20)] ...
                            = deal(NaN);

%% EVALUATE JORDAN/BOWLING/LEYSINGER-VIELI/PANTON+KARLSSON DATA    

disp('Generate basal water mask from Jordan et al. (2018), Bowling et al. (2019), Panton and Karlsson (2015) and Leysinger-Vieli et al. (2018) data...')

% mask initializations
[mask_basal_water_J18, mask_basal_water_B19, mask_plume_LV18, mask_plume_PK15] ...
                            = deal(zeros(num_decim_y, num_decim_x));
[mask_basal_water_J18(~mask_gris_decim), mask_basal_water_B19(~mask_gris_decim), mask_plume_LV18(~mask_gris_decim), mask_plume_PK15(~mask_gris_decim)] ...
                            = deal(NaN);

% Jordan et al. (2018) basal water identifications
ind_J18_xy_near             = [interp1(x_decim(1, :), 1:num_decim_x, (1e-3 .* J18.x(find(J18.basal_water))), 'nearest', 'extrap') interp1(y_decim(:, 1)', 1:num_decim_y, (1e-3 .* J18.y(find(J18.basal_water))), 'nearest', 'extrap')];
ind_J18_xy_unique           = unique(ind_J18_xy_near, 'rows'); % unique cells with basal water IDs
% given numeric value to number of times water found in a cell
for ii = 1:size(ind_J18_xy_unique, 1)
    mask_basal_water_J18(ind_J18_xy_unique(ii, 2), ind_J18_xy_unique(ii, 1)) ...
                            = length(find(ismember(ind_J18_xy_near, ind_J18_xy_unique(ii, :), 'rows')));
end

% Bowling et al. (2019) subglacial lakes
ind_B19_norepeat            = find(~isnan(B19.Lat)); % don't need repeat lakes
ind_B19_xy_near             = [interp1(x_decim(1, :), 1:num_decim_x, (1e-3 .* B19.x(ind_B19_norepeat)), 'nearest', 'extrap') interp1(y_decim(:, 1)', 1:num_decim_y, (1e-3 .* B19.y(ind_B19_norepeat)), 'nearest', 'extrap')];

% assign numeric confidence levels to Bowling et al. (2019) confidence assessments
B19.ConfidenceLevelNum      = zeros(size(B19, 1), 1);
B19.ConfidenceLevelNum(B19.ConfidenceLevel == 'Low') ...
                            = 1;
B19.ConfidenceLevelNum(B19.ConfidenceLevel == 'Medium') ...
                            = 5;
B19.ConfidenceLevelNum(B19.ConfidenceLevel == 'High') ...
                            = 9;
B19.ConfidenceLevelNum(B19.ConfidenceLevel == 'Very high') ...
                            = 10;
B19.ConfidenceLevelNum(intersect(find(B19.ConfidenceLevelNum == 0), ind_B19_norepeat)) ...
                            = 5; % undefined confidence level gets medium
for ii = 1:size(ind_B19_xy_near, 1)
    mask_basal_water_B19(ind_B19_xy_near(ii, 2), ind_B19_xy_near(ii, 1)) ...
                            = mask_basal_water_B19(ind_B19_xy_near(ii, 2), ind_B19_xy_near(ii, 1)) + B19.ConfidenceLevelNum(ind_B19_norepeat(ii)); % add confidence value to cell
end

% Leysinger-Vieli et al. (2018) basal plume identifications
% start with small plumes to which we assign unity
ind_LV18_small_xy_near      = [interp1(x_decim(1, :), 1:num_decim_x, (1e-3 .* LV18.x_plume_s), 'nearest', 'extrap') interp1(y_decim(:, 1)', 1:num_decim_y, (1e-3 .* LV18.y_plume_s), 'nearest', 'extrap')];
ind_LV18_small_xy_unique    = unique(ind_LV18_small_xy_near, 'rows');
for ii = 1:size(ind_LV18_small_xy_unique, 1)
    mask_plume_LV18(ind_LV18_small_xy_unique(ii, 2), ind_LV18_small_xy_unique(ii, 1)) ...
                            = length(find(ismember(ind_LV18_small_xy_near, ind_LV18_small_xy_unique(ii, :), 'rows')));
end
% add in large plumes (value=5)
ind_LV18_large_xy_near      = [interp1(x_decim(1, :), 1:num_decim_x, (1e-3 .* LV18.x_plume_l), 'nearest', 'extrap') interp1(y_decim(:, 1)', 1:num_decim_y, (1e-3 .* LV18.y_plume_l), 'nearest', 'extrap')];
ind_LV18_large_xy_unique    = unique(ind_LV18_large_xy_near, 'rows');
for ii = 1:size(ind_LV18_large_xy_unique, 1)
    mask_plume_LV18(ind_LV18_large_xy_unique(ii, 2), ind_LV18_large_xy_unique(ii, 1)) ...
                            = mask_plume_LV18(ind_LV18_large_xy_unique(ii, 2), ind_LV18_large_xy_unique(ii, 1)) + (5 * length(find(ismember(ind_LV18_small_xy_near, ind_LV18_small_xy_unique(ii, :), 'rows'))));
end

% Panton and Karlsson (2015) units of disrupted radiostratigraphy (UDRs)
thick_UDR                   = interp2(BM4.x, BM4.y, BM4.thick, x_UDR, y_UDR); % BM4 ice thickess at UDRs
thick_UDR_rel               = PK15.UDR ./ thick_UDR; % relative thickness of UDRs
thick_UDR_min               = 1e3; % minimum absolute ice thickness for which to consider UDR
thick_UDR_rel_cutoff        = (1 / 3); % cutoff above which to assume a "large" plume, following methods of Leysinger-Vieli et al. (2018)
% start with small UDRs to which we assign unity
ind_PK15_small_xy_near      = [interp1(x_decim(1, :), 1:num_decim_x, x_UDR((thick_UDR_rel < thick_UDR_rel_cutoff) & (thick_UDR >= thick_UDR_min)), 'nearest', 'extrap') ...
                               interp1(y_decim(:, 1)', 1:num_decim_y, y_UDR((thick_UDR_rel < thick_UDR_rel_cutoff) & (thick_UDR >= thick_UDR_min)), 'nearest', 'extrap')];
ind_PK15_small_xy_unique    = unique(ind_PK15_small_xy_near, 'rows');
for ii = 1:size(ind_PK15_small_xy_unique, 1)
    mask_plume_PK15(ind_PK15_small_xy_unique(ii, 2), ind_PK15_small_xy_unique(ii, 1)) ...
                            = length(find(ismember(ind_PK15_small_xy_near, ind_PK15_small_xy_unique(ii, :), 'rows')));
end
% add in large UDRs (value=5)
ind_PK15_large_xy_near      = [interp1(x_decim(1, :), 1:num_decim_x, x_UDR((thick_UDR_rel >= thick_UDR_rel_cutoff) & (thick_UDR >= thick_UDR_min)), 'nearest', 'extrap') ...
                               interp1(y_decim(:, 1)', 1:num_decim_y, y_UDR((thick_UDR_rel >= thick_UDR_rel_cutoff) & (thick_UDR >= thick_UDR_min)), 'nearest', 'extrap')];
ind_PK15_large_xy_unique    = unique(ind_PK15_large_xy_near, 'rows');
for ii = 1:size(ind_PK15_large_xy_unique, 1)
    mask_plume_PK15(ind_PK15_large_xy_unique(ii, 2), ind_PK15_large_xy_unique(ii, 1)) ...
                            = mask_plume_PK15(ind_PK15_large_xy_unique(ii, 2), ind_PK15_large_xy_unique(ii, 1)) + (5 * length(find(ismember(ind_PK15_small_xy_near, ind_PK15_small_xy_unique(ii, :), 'rows'))));
end

% combine plume masks from both studies (PK15 and LV18)
mask_plume_combo            = max(cat(3, mask_plume_LV18, mask_plume_PK15), [], 3);

% total mask including all basal water identifications by sum
mask_basal_water            = mask_basal_water_J18 + mask_basal_water_B19 + mask_plume_combo;

% expand each mask by 1 cell in each direction (including diagonal)
tmp                         = mask_basal_water_J18; % start with initial mask
tmp(isnan(tmp))             = 0; % get rid of NaNs
mask_basal_water_J18_exp    = imdilate(tmp, strel('square', 3)); % _exp for expansion
mask_basal_water_J18_exp(~mask_gris_decim) ...
                            = NaN; % re-NaN
tmp                         = mask_basal_water_B19;
tmp(isnan(tmp))             = 0;
mask_basal_water_B19_exp    = imdilate(tmp, strel('square', 3));
mask_basal_water_B19_exp(~mask_gris_decim) ...
                            = NaN;
tmp                         = mask_plume_combo;
tmp(isnan(tmp))             = 0;
mask_plume_combo_exp        = imdilate(tmp, strel('square', 3));
mask_plume_combo_exp(~mask_gris_decim) ...
                            = NaN;
tmp                         = mask_basal_water;
tmp(isnan(tmp))             = 0;
mask_basal_water            = imdilate(tmp, strel('square', 3));
mask_basal_water(~mask_gris_decim) ...
                            = NaN;

% thresholds for later agreememt calculations
threshold_basal_water_std   = 5;
threshold_basal_water_cold  = 10;
threshold_basal_water_warm  = 1;

% build simpler masks for agreement calculation
[mask_basal_water_std, mask_basal_water_cold, mask_basal_water_warm] ...
                            = deal(zeros(num_decim_y, num_decim_x)); % initialize to 0
[mask_basal_water_std(~mask_gris_decim), mask_basal_water_cold(~mask_gris_decim), mask_basal_water_warm(~mask_gris_decim)] ...
                            = deal(NaN); % NaN out non-ice
mask_basal_water_std(mask_basal_water >= threshold_basal_water_std) ...
                            = 1; % unity / logical for positive agreement
mask_basal_water_cold(mask_basal_water >= threshold_basal_water_cold) ...
                            = 1;
mask_basal_water_warm(mask_basal_water >= threshold_basal_water_warm) ...
                            = 1;

%% ISMIP6

disp('Loading ISMIP6 data and assessing agreement...')

% load ISMIP6 .mat compilation of all ctrl/ctrl_proj runs
ismip6                      = load('mat/ismip6');
[ismip6.temp_bed_interp, ismip6.temp_bed_pmp, ismip6.thick_interp] ...
                            = deal(cell(1, ismip6.num_inst));
for ii = 1:ismip6.num_inst
    [ismip6.temp_bed_interp{ii}, ismip6.temp_bed_pmp{ii}, ismip6.thick_interp{ii}] ...
                            = deal(cell(1, ismip6.num_model(ii)));    
    for jj = 1:ismip6.num_model(ii)
        if ((ii == 5) && jj == 2)
            ismip6.thick{ii}{jj}(ismip6.thick{ii}{jj} <= 2) ...
                            = NaN; % get rid of silly values
        end        
        ismip6.thick_interp{ii}{jj} ...
                            = interp2(ismip6.x, ismip6.y, double(ismip6.thick{ii}{jj}), x_decim, y_decim, 'nearest'); % model thickness at decimated grid
        ismip6.temp_bed_interp{ii}{jj} ...
                            = interp2(ismip6.x, ismip6.y, ismip6.temp_bed{ii}{jj}, x_decim, y_decim, 'nearest'); % nearest-neighbor interpolated basal temperature
        ismip6.temp_bed_pmp{ii}{jj} ...
                            = ismip6.temp_bed_interp{ii}{jj} + (273.15 - pmp(ismip6.thick_interp{ii}{jj})); % pressure melting correction
        ismip6.temp_bed_pmp{ii}{jj}(ismip6.temp_bed_pmp{ii}{jj} > 0) ...
                            = 0; % fix over-corrections
    end
end

% downselect ISMIP6 model instances to include
ismip6_incl_agree           = [1  3; % AWI/ISSM3
                               4  2; % ILTS_PIK/SICOPOLIS3
                               5  2; % JPL/ISSMPALEO
                               6  2; % LSCE/GRISLI2
                               7  1; % MUN/GSM2601
                               8  1; % NCAR/CISM
                               9  2; % UAF/PISM2
                               10 1; % UCIJPL/ISSM1
                               11 1; % VUW/GISMHOMv1
                               12 1]; % VUW/PISM
num_incl_agree              = size(ismip6_incl_agree, 1);

% initialize ISMIP6 agreement masks
[mask_agree_ismip6, mask_agree_ismip6_cold, mask_agree_ismip6_warm, num_model_mat_ismip6] ...
                            = deal(NaN(num_decim_y, num_decim_x));
[mask_agree_ismip6(mask_gris_decim), mask_agree_ismip6_cold(mask_gris_decim), mask_agree_ismip6_warm(mask_gris_decim), num_model_mat_ismip6(mask_gris_decim)] ...
                            = deal(0);

% temperature thresholds for thawed, effectively confidence level in prediction, degrees Celsius
temp_pmp_std                = -1.0;
temp_pmp_cold               = -0.5;
temp_pmp_warm               = -1.5;

% ISMIP6 agreement masks
for ii = 1:num_incl_agree
    [jj, kk]                = deal(ismip6_incl_agree(ii, 1), ismip6_incl_agree(ii, 2));
    num_model_mat_ismip6(~isnan(ismip6.temp_bed_pmp{jj}{kk})) ...
                            = num_model_mat_ismip6(~isnan(ismip6.temp_bed_pmp{jj}{kk})) + 1;
    mask_agree_ismip6(ismip6.temp_bed_pmp{jj}{kk} >= temp_pmp_std) ...
                            = mask_agree_ismip6(ismip6.temp_bed_pmp{jj}{kk} >= temp_pmp_std) + 1;
    mask_agree_ismip6(ismip6.temp_bed_pmp{jj}{kk} < temp_pmp_std) ...
                            = mask_agree_ismip6(ismip6.temp_bed_pmp{jj}{kk} < temp_pmp_std) - 1;
    mask_agree_ismip6_cold(ismip6.temp_bed_pmp{jj}{kk} >= temp_pmp_cold) ...
                            = mask_agree_ismip6_cold(ismip6.temp_bed_pmp{jj}{kk} >= temp_pmp_cold) + 1;
    mask_agree_ismip6_cold(ismip6.temp_bed_pmp{jj}{kk} < temp_pmp_cold) ...
                            = mask_agree_ismip6_cold(ismip6.temp_bed_pmp{jj}{kk} < temp_pmp_cold) - 1;
    mask_agree_ismip6_warm(ismip6.temp_bed_pmp{jj}{kk} >= temp_pmp_warm) ...
                            = mask_agree_ismip6_warm(ismip6.temp_bed_pmp{jj}{kk} >= temp_pmp_warm) + 1;
    mask_agree_ismip6_warm(ismip6.temp_bed_pmp{jj}{kk} < temp_pmp_warm) ...
                            = mask_agree_ismip6_warm(ismip6.temp_bed_pmp{jj}{kk} < temp_pmp_warm) - 1;
end

% normalize mask but number available at each location
[mask_agree_ismip6, mask_agree_ismip6_cold, mask_agree_ismip6_warm] ...
                            = deal((mask_agree_ismip6 ./ num_incl_agree), (mask_agree_ismip6_cold ./ num_incl_agree), (mask_agree_ismip6_warm ./ num_incl_agree));

mask_agree_diff             = mask_agree_ismip6 - GBaTSv1.mask_agree_searise_update; % difference between ISMIP6 and SeaRISE basal thermal state agreement masks

% bed elevation at decimated mask needed for ISMIP6 display, better to use decimated not filtered
BM4.elev_bed_decim         = interp2(BM4.x, BM4.y, BM4.elev_bed, x_decim, y_decim);

%% GBaTSv2 AGREEMENT AND LIKELY STATE MASKS

disp('Calculate agreement and likely basal thermal state masks...')

[mask_agree, mask_agree_cold, mask_agree_warm, mask_likely, mask_num_constraint, mask_num_constraint_cold, mask_num_constraint_warm] ...
                            = deal(NaN(num_decim_y, num_decim_x));
[mask_agree(mask_gris_decim), mask_agree_cold(mask_gris_decim), mask_agree_warm(mask_gris_decim), mask_likely(mask_gris_decim), mask_num_constraint(mask_gris_decim), mask_num_constraint_cold(mask_gris_decim), mask_num_constraint_warm(mask_gris_decim)] ...
                            = deal(0);

% basal thermal state constraints
name_constraint             = {'melt_bed_filt'       'speed_ratio_std_trim'      'mask_agree_ismip6'          'mask_basal_water_std'};
name_constraint_cold        = {'melt_bed_lo_filt'    'speed_ratio_min_trim'      'mask_agree_ismip6_cold'     'mask_basal_water_cold'};
name_constraint_warm        = {'melt_bed_hi_filt'    'speed_ratio_max_trim'      'mask_agree_ismip6_warm'     'mask_basal_water_warm'};

% numeric constraints
constraint                  = [NaN 0.01;  NaN 1;              -0.3 0.3;                    NaN 1];
num_constraint              = length(name_constraint);

num_hole_small              = 10; % number of pixels for hole threshold

% loop through each constraint, evaluate hole-filled agreement, +1 if thawed, -1 if frozen
for ii = 1:num_constraint
    ind_std                 = find(smallholefill(eval(name_constraint{ii}), num_hole_small, constraint(ii, 2), '>=') >= constraint(ii, 2)); % points in hole-filled mask that meet threshold
    [mask_agree(ind_std), mask_num_constraint(ind_std)] ...
                            = deal((mask_agree(ind_std) + 1), (mask_num_constraint(ind_std) + 1));
    ind_cold                = find(smallholefill(eval(name_constraint_cold{ii}), num_hole_small, constraint(ii, 2), '>=') >= constraint(ii, 2));    
    [mask_agree_cold(ind_cold), mask_num_constraint_cold(ind_cold)] ...
                            = deal((mask_agree_cold(ind_cold) + 1), (mask_num_constraint_cold(ind_cold) + 1));
    ind_warm                = find(smallholefill(eval(name_constraint_warm{ii}), num_hole_small, constraint(ii, 2), '>=') >= constraint(ii, 2));    
    [mask_agree_warm(ind_warm), mask_num_constraint_warm(ind_warm)] ...
                            = deal((mask_agree_warm(ind_warm) + 1), (mask_num_constraint_warm(ind_warm) + 1));
    if ~isnan(constraint(ii, 1)) % cold constraint available, only available for ISMIP6 comparison
        ind_std             = find(smallholefill(eval(name_constraint{ii}), num_hole_small, constraint(ii, 1), '<=') <= constraint(ii, 1));        
        [mask_agree(ind_std), mask_num_constraint(ind_std)] ...
                            = deal((mask_agree(ind_std) - 1), (mask_num_constraint(ind_std) + 1));
        ind_cold            = find(smallholefill(eval(name_constraint_cold{ii}), num_hole_small, constraint(ii, 1), '<=') <= constraint(ii, 1));        
        [mask_agree_cold(ind_cold), mask_num_constraint_cold(ind_cold)] ...
                            = deal((mask_agree_cold(ind_cold) - 1), (mask_num_constraint_cold(ind_cold) + 1));
        ind_warm            = find(smallholefill(eval(name_constraint_warm{ii}), num_hole_small, constraint(ii, 1), '<=') <= constraint(ii, 1));        
        [mask_agree_warm(ind_warm), mask_num_constraint_warm(ind_warm)] ...
                            = deal((mask_agree_warm(ind_warm) - 1), (mask_num_constraint_warm(ind_warm) + 1));
    end
end

% normalize synthesis agreement masks to within +/- 1
[mask_agree, mask_agree_cold, mask_agree_warm] ...
                            = deal((mask_agree ./ num_constraint), (mask_agree_cold ./ num_constraint), (mask_agree_warm ./ num_constraint));

% likely basal thermal state mask (GBaTSv2)
mask_likely((mask_agree > 0) & (mask_agree_cold >= 0) & (mask_agree_warm > 0)) ...
                            = 1; % likely thawed
mask_likely((mask_agree < 0) & (mask_agree_cold < 0) & (mask_agree_warm <= 0)) ...
                            = -1; % likely frozen

mask_likely_filled          = smallholefill(smallholefill(mask_likely, num_hole_small, 1, '=='), num_hole_small, -1, '=='); % fill likely thawed then frozen small holes

% statistics 
numel_good                  = length(find(mask_gris_decim & (thick_filt > 0))); % number of on-ice thick elements

% ISMIP6 frozen
frac_ismip6_frozen          = round(1e2 * (length(find(mask_agree_ismip6 <= constraint(3, 1))) / numel_good));
frac_ismip6_frozen_cold     = round(1e2 * (length(find(mask_agree_ismip6_cold <= constraint(3, 1))) / numel_good));
frac_ismip6_frozen_warm     = round(1e2 * (length(find(mask_agree_ismip6_warm <= constraint(3, 1))) / numel_good));

% ISMIP6 thawed
frac_ismip6_thawed          = round(1e2 * (length(find(mask_agree_ismip6 >= constraint(3, 2))) / numel_good));
frac_ismip6_thawed_cold     = round(1e2 * (length(find(mask_agree_ismip6_cold >= constraint(3, 2))) / numel_good));
frac_ismip6_thawed_warm     = round(1e2 * (length(find(mask_agree_ismip6_warm >= constraint(3, 2))) / numel_good));

% ISMIP6 uncertainty
frac_ismip6_uncert          = round(1e2 * (length(find(abs(mask_agree_ismip6) < constraint(3, 2))) / numel_good));
frac_ismip6_uncert_cold     = round(1e2 * (length(find(abs(mask_agree_ismip6_cold) < constraint(3, 2))) / numel_good));
frac_ismip6_uncert_warm     = round(1e2 * (length(find(abs(mask_agree_ismip6_warm) < constraint(3, 2))) / numel_good));

% surface/deformation speed thawed
frac_speed_thawed           = round(1e2 * (length(find(speed_ratio_std_trim >= constraint(2, 2))) / numel_good));
frac_speed_thawed_min       = round(1e2 * (length(find(speed_ratio_min_trim >= constraint(2, 2))) / numel_good));
frac_speed_thawed_max       = round(1e2 * (length(find(speed_ratio_max_trim >= constraint(2, 2))) / numel_good));

% basal melt frozen/thawed
frac_melt_layer_thawed      = round(1e2 * (length(find(melt_bed_filt(mask_D1_decim) >= constraint(1, 2))) / numel_good));
frac_melt_layer_thawed_min  = round(1e2 * (length(find(melt_bed_lo_filt(mask_D1_decim) >= constraint(1, 2))) / numel_good));
frac_melt_layer_thawed_max  = round(1e2 * (length(find(melt_bed_hi_filt(mask_D1_decim) >= constraint(1, 2))) / numel_good));

% frozen/thawed
frac_melt_id_thawed         = round(1e2 * (length(find(mask_basal_water_std >= constraint(1, 2))) / numel_good));
frac_melt_id_thawed_min     = round(1e2 * (length(find(mask_basal_water_cold >= constraint(1, 2))) / numel_good));
frac_melt_id_thawed_max     = round(1e2 * (length(find(mask_basal_water_warm >= constraint(1, 2))) / numel_good));

% likely basal thermal state mask
frac_frozen                 = round(1e2 * (length(find(mask_likely == -1)) / numel_good));
frac_thawed                 = round(1e2 * (length(find(mask_likely == 1)) / numel_good));
frac_uncert                 = round(1e2 * (length(find(mask_likely == 0)) / numel_good));

% likely mask basal thermal state mask, hole-filled
frac_frozen_filled          = round(1e2 * (length(find(mask_likely_filled == -1)) / numel_good));
frac_thawed_filled          = round(1e2 * (length(find(mask_likely_filled == 1)) / numel_good));
frac_uncert_filled          = round(1e2 * (length(find(mask_likely_filled == 0)) / numel_good));

%% basin areas and BTS fractions

disp('Assign basin areas to 5-km grid and compare v1/v2...')

ind_basin                   = NaN(num_decim_y, num_decim_x);
for ii = 1:num_M19
    ind_basin_tmp           = inpolygon(x_decim, y_decim, M19(ii).X, M19(ii).Y);
    ind_basin(ind_basin_tmp)= ii;
end

% basal thermal state by basin for both v1 and v2
[basin_area_GBaTSv1, basin_area_GBaTSv2] ...
                            = deal(NaN(num_M19, 3));
for ii = 1:num_M19
    basin_area_GBaTSv1(ii, :) ...
                            = round(1e2 * ([length(find((GBaTSv1.mask_likely_filled == -1) & (ind_basin == ii))) length(find((GBaTSv1.mask_likely_filled == 0) & (ind_basin == ii))) length(find((GBaTSv1.mask_likely_filled == 1) & (ind_basin == ii)))] ./ ...
                                           length(find((ind_basin == ii) & (thick_filt > 0)))));    
    basin_area_GBaTSv2(ii, :) ...
                            = round(1e2 * ([length(find((mask_likely_filled == -1) & (ind_basin == ii))) length(find((mask_likely_filled == 0) & (ind_basin == ii))) length(find((mask_likely_filled == 1) & (ind_basin == ii)))] ./ length(find((ind_basin == ii) & (thick_filt > 0)))));
end

%% SeaRISE imports

disp('Load SeaRISE ice thickness used for comparison against ISMIP6...')

SeaRISE                     = struct;
SeaRISE.x                   = double(ncread('/Users/jamacgre/Documents/data/GBaTSv2/SeaRISE/Greenland_5km_E2.nc', 'x1'))';
SeaRISE.y                   = double(ncread('/Users/jamacgre/Documents/data/GBaTSv2/SeaRISE/Greenland_5km_E2.nc', 'y1'));
SeaRISE.elev_bed_orig       = double(ncread('/Users/jamacgre/Documents/data/GBaTSv2/SeaRISE/Greenland_5km_E2.nc', 'topg'))';

wgs84                       = wgs84Ellipsoid;

% convert to 5 km grid
[SeaRISE.x_grd, SeaRISE.y_grd] ...
                            = meshgrid(SeaRISE.x, SeaRISE.y);
[SeaRISE.lat, SeaRISE.lon]  = polarstereo_inv(SeaRISE.x_grd, SeaRISE.y_grd, wgs84.SemimajorAxis, wgs84.Eccentricity, 71, -39);

[SeaRISE.x_tmp, SeaRISE.y_tmp] ...
                            = projfwd(projcrs(3413), SeaRISE.lat, SeaRISE.lon);
[SeaRISE.x_tmp, SeaRISE.y_tmp] ...
                            = deal((1e-3 .* SeaRISE.x_tmp(:)), (1e-3 .* SeaRISE.y_tmp(:)));
SeaRISE.elev_bed_interp     = scatteredInterpolant(SeaRISE.x_tmp, SeaRISE.y_tmp, SeaRISE.elev_bed_orig(:), 'natural', 'none');
SeaRISE.elev_bed            = SeaRISE.elev_bed_interp(x_decim, y_decim);
SeaRISE.thick               = (interp2(BM3.x, BM3.y, BM3.elev_surf, x_decim, y_decim) - SeaRISE.elev_bed); % interpolated SeaRISE thickness onto EPSG:3413 5 km grid

BM3.thick_interp            = interp2(BM3.x, BM3.y, BM3.thick, x_decim, y_decim); % BM3 thickness on 5-km grid

thick_change                = BM3.thick_interp - SeaRISE.thick; % thickness change from SeaRISE to ISMIP6
ind_good                    = find(mask_gris_decim & ~isnan(thick_change) & ~isnan(mask_agree_diff)); % indices where differences relevant/sane
model_change_bin            = hist3([thick_change(ind_good) (1e2 .* mask_agree_diff(ind_good))], 'Edges', {-2000:100:2000 -150:10:150})' .* (1e2 / length(ind_good)); % binned changes in thickness and basal temperature agreement

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
    
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
   
%% FIGURE 1: BOREHOLES / LAKES OVERLAIN ON DRIVING STRESS
    
    figure('position', [200 200 500 800], 'color', 'w', 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; parula(10)])
    subplot('position', [0.04 0.03 0.9 0.92]);
    hold on
    range_tmp               = linspace(0, 150, 11);
    plot_tmp                = interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 4:13, (1e-3 .* driving_stress), 'nearest', 'extrap');
    plot_tmp(isnan(driving_stress)) ...
                            = mask_combo_decim(isnan(driving_stress)) + 1;
    imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
    for ii = 1:num_coast
        plot(x_coast{ii}, y_coast{ii}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
    end
    for ii = 1:length(x_paral)
        plot((1e-3 .* x_paral{ii}), (1e-3 .* y_paral{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:length(x_merid)
        plot((1e-3 .* x_merid{ii}), (1e-3 .* y_merid{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:num_M19
        plot(M19(ii).X, M19(ii).Y, 'k', 'linewidth', 2.5)
    end
    pr                      = plot(NaN, NaN, 'k', 'linewidth', 1);
    pbf                     = plot(x_borehole_frozen, y_borehole_frozen, 'wo', 'markersize', 16, 'markerfacecolor', [0 0 0.75], 'linewidth', 1);
    pbt                     = plot(x_borehole_thawed, y_borehole_thawed, 'wo', 'markersize', 16, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
    psl                     = plot(x_lake, y_lake, 'wd', 'markersize', 14, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
    tb                      = NaN(1, num_M19);
    tb_pos                  = 1e3 .* [-0.0011   -2.1028         0
                                       0.0964   -2.5597         0
                                      -0.1403   -2.5748         0
                                      -0.2850   -1.5266         0
                                      -0.0751   -1.1437         0
                                       0.2128   -1.4003         0
                                       0.3049   -2.0452         0];
    for ii = 1:num_M19
        tb(ii)              = text(mean(M19(ii).X), mean(M19(ii).Y), M19(ii).SUBREGION1, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold');
        set(tb(ii), 'position', tb_pos(ii, :))
    end
    fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
    fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
    text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
    text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
    text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
    text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
    text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
    text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
    caxis([1 14])
    axis equal
    axis([x_min x_max y_min y_max])
    box on
    grid on
    set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
    lg                      = legend([pbf pbt psl], 'borehole (frozen bed)', 'borehole (thawed bed)', 'subglacial lake');
    set(lg, 'location', 'northoutside', 'color', [0.9 0.9 0.9])
    tick_str                = cell(1, 21);
    for ii = 1:11
        tick_str{ii}        = num2str(range_tmp(ii));
    end
    tick_str{end}           = ['>' tick_str{end}];
    colorbar('ylim', [4 14], 'ytick', 4:14, 'yticklabel', tick_str, 'fontsize', 20, 'ticklength', 0.0225, 'fontweight', 'bold')
    text(885, -470, '$\tau_d$', 'fontsize', 20, 'interpreter', 'latex', 'fontweight', 'bold')
    text(850, -570, '(kPa)', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold')

%% FIGURE 2: IMSIP6 MODEL COMPARISON (BASAL TEMPERATURE AT END OF EACH CONTROL RUN)
   
    load('/Users/jamacgre/Documents/research/matlab/greenland/mat/rb10', 'rb10')
    cmap1                   = [([135 206 235] ./ 255); ([159 89 39] ./ 255); 1 1 1];
    cmap2                   = rb10;
    figure('position', [100 200 1080 720], 'color', 'w', 'renderer', 'zbuffer')
    colormap([cmap1; cmap2])
    ax                      = NaN(2, 5);
    ranges                  = [-10 0];
    incs                    = 1;
    range_tmp               = ranges(1):incs:ranges(2);
    tick_str                = cell(1, length(range_tmp));
    for ii = 1:length(range_tmp)
        tick_str{ii}        = num2str(range_tmp(ii));
    end
    tick_str{1}             = ['< ' tick_str{1}];
    ind_ord                 = [1 3 8 7 10 5 9 2 4 6];
    for ii = ind_ord
        if (ii <= 5)
            axes('position', [(((ii - 1) * (0.32 / ((y_max - y_min) / (x_max - x_min)))) + 0.01) 0.49 (0.30 / ((y_max - y_min) / (x_max - x_min))) 0.50]);
        else
            axes('position', [(((ii - 6) * (0.32 / ((y_max - y_min) / (x_max - x_min)))) + 0.01) 0.00 (0.30 / ((y_max - y_min) / (x_max - x_min))) 0.50]);
        end
        hold on
        plot_tmp            = 3 + discretize(ismip6.temp_bed_pmp{ismip6_incl_agree(ind_ord(ii), 1)}{ismip6_incl_agree(ind_ord(ii), 2)}, [-Inf range_tmp(2:(end - 1)) Inf]);
        plot_tmp(isnan(plot_tmp) & ((mask_combo_decim == 0) | ((mask_combo_decim == 2) & (BM4.elev_bed_decim  <= 0)))) ...
                            = 1;
        plot_tmp(isnan(plot_tmp) & ((mask_combo_decim == 1) | ((mask_combo_decim == 2) & (BM4.elev_bed_decim  > 0)))) ...
                            = 2;
        imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 0.5)
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 0.5)
        end
        for jj = 1:length(x_borehole_frozen)
            plot(x_borehole_frozen(jj), y_borehole_frozen(jj), 'wo', 'markersize', 8, 'markerfacecolor', cmap2(discretize(temp_borehole_frozen(jj), [-Inf range_tmp(2:(end - 1)) Inf]), :), 'linewidth', 1);
        end        
        plot(x_borehole_thawed, y_borehole_thawed, 'wo', 'markersize', 8, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
        plot(x_borehole_thawed(10), y_borehole_thawed(10), 'wo', 'markersize', 8, 'markerfacecolor', cmap2(discretize(-1.3, [-Inf range_tmp(2:(end - 1)) Inf]), :), 'linewidth', 1);
        plot(x_lake, y_lake, 'wd', 'markersize', 7, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
        if (ii == 1)
            fill([325 575 575 325], [-3200 -3200 -3250 -3250], 'k')
            text(315, -3140, '250 km', 'color', 'k', 'fontsize', 10, 'fontweight', 'bold')
            text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 10, 'rotation', -10)
            text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 10, 'rotation', 82)
            text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 10, 'rotation', 12)
            text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 10, 'rotation', -75)
        end
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
        title(['(' letters(ii) ') ' ismip6.name_inst{ismip6_incl_agree(ind_ord(ii), 1)} '/' ismip6.name_model{ismip6_incl_agree(ind_ord(ii), 1)}{ismip6_incl_agree(ind_ord(ii), 2)}], 'color', 'k', 'fontweight', 'bold', 'fontsize', 16, 'interpreter', 'none')
        caxis([1 (3 + length(range_tmp))])
        if (ii == num_incl_agree)
            colorbar('position', [0.90 0.025 0.015 0.94], 'fontsize', 16, 'ylim', [4 (3 + length(range_tmp))], 'ytick', 4:(3 + length(range_tmp)), 'yticklabel', tick_str, 'ticklength', 0.025, 'fontweight', 'bold', 'ticklength', 0.0225, 'fontsize', 16)
            text(1433, 1537, 'Pressure-corrected basal temperature at end of control run (\circC)', 'color', 'k', 'fontsize', 16, 'fontweight', 'bold', 'rotation', 270)
        end
    end
    
%% FIGURE 3: ISMIP6 AGREEMENT MAP
    
    figure('position', [200 200 (1.5 * 500) (((y_max - y_min) / (x_max - x_min)) * 500)], 'color', 'w', 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; redblue(11, 0.5)])
    subplot('position', [-0.01 0.02 0.8 0.95])
    hold on
    range_tmp               = linspace(-1, 1, (num_incl_agree + 1));
    plot_tmp                = NaN(num_decim_y, num_decim_x);
    plot_tmp(mask_agree_ismip6 < 0) ...
                            = discretize(mask_agree_ismip6(mask_agree_ismip6 < 0), range_tmp, 'includededge', 'left');
    plot_tmp(mask_agree_ismip6 == 0) ...
                            = find(range_tmp == 0);
    plot_tmp(mask_agree_ismip6 > 0) ...
                            = discretize(mask_agree_ismip6(mask_agree_ismip6 > 0), range_tmp, 'includededge', 'right');
    plot_tmp                = 4 + plot_tmp;
    plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
    imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
    for ii = 1:num_coast
        plot(x_coast{ii}, y_coast{ii}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
    end
    for ii = 1:length(x_paral)
        plot((1e-3 .* x_paral{ii}), (1e-3 .* y_paral{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:length(x_merid)
        plot((1e-3 .* x_merid{ii}), (1e-3 .* y_merid{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:num_M19
        plot(M19(ii).X, M19(ii).Y, 'k', 'linewidth', 2)
    end
    pbf                     = plot(x_borehole_frozen, y_borehole_frozen, 'wo', 'markersize', 16, 'markerfacecolor', [0 0 0.75], 'linewidth', 1);
    pbt                     = plot(x_borehole_thawed, y_borehole_thawed, 'wo', 'markersize', 16, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
    psl                     = plot(x_lake, y_lake, 'wd', 'markersize', 14, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
    fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
    fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
    text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
    text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
    text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
    text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
    text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 13)
    text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
    caxis([1 (5 + num_incl_agree)])
    axis equal
    axis([x_min x_max y_min y_max])
    box on
    set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
    colorbar('fontsize', 20, 'location', 'eastoutside', 'limits', [4 (5 + num_incl_agree)], 'ytick', 4:(5 + num_incl_agree), 'ticklength', 0.018, 'fontweight', 'bold', 'yticklabel', {})
    tick_str                = {'all agree frozen bed' '9 / 10 frozen' '8 / 10 frozen' '7 / 10 frozen' '6 / 10 frozen' '5 frozen, 5 thawed' '6 / 10 thawed' '7 / 10 thawed' '8 / 10 thawed' '9 / 10 thawed' 'all agree thawed bed'};
    for ii = 1:(num_incl_agree + 1)
        text(940, (-3475 + (ii * 245)), tick_str{ii}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold')
    end   
    
%% FIGURE 4: NYE+MELT MELT RATE
    
    rb                      = redblue(20, 0.5);
    figure('position', [200 200 (1.25 * 500) (((y_max - y_min) / (x_max - x_min)) * 500)], 'color', 'w', 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; rb([9 11:end], :)])
    subplot('position', [0.05 0.02 0.8 0.90])
    hold on
    range_tmp               = linspace(-0.01, 0.1, 12);
    plot_tmp                = 3 + discretize(melt_bed_filt, [-Inf range_tmp(2:(end - 1)) Inf]);
    plot_tmp(isnan(plot_tmp) | ~mask_D1_decim) ...
                            = mask_combo_decim(isnan(plot_tmp) | ~mask_D1_decim) + 1;
    imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
    for ii = 1:num_coast
        plot(x_coast{ii}, y_coast{ii}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
    end
    for ii = 1:length(x_paral)
        plot((1e-3 .* x_paral{ii}), (1e-3 .* y_paral{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:length(x_merid)
        plot((1e-3 .* x_merid{ii}), (1e-3 .* y_merid{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:num_M19
        plot(M19(ii).X, M19(ii).Y, 'k', 'linewidth', 2)
    end
    pbf                     = plot(x_borehole_frozen, y_borehole_frozen, 'wo', 'markersize', 16, 'markerfacecolor', [0 0 0.75], 'linewidth', 1);
    pbt                     = plot(x_borehole_thawed, y_borehole_thawed, 'wo', 'markersize', 16, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
    psl                     = plot(x_lake, y_lake, 'wd', 'markersize', 14, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
    fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
    fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
    text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
    text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
    text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
    text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
    text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 13)
    text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
    caxis([1 15])
    axis equal
    axis([x_min x_max y_min y_max])
    box on
    set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
    tick_str                = cell(1, 11);
    for ii = 1:12
        tick_str{ii}        = num2str(1e2 * range_tmp(ii));
    end
    tick_str{1}             = ['<' tick_str{1}];
    tick_str{end}           = ['>' tick_str{end}];    
    colorbar('fontsize', 20, 'location', 'eastoutside', 'limits', [4 15], 'ytick', 4:15, 'yticklabel', tick_str, 'ticklength', 0.018, 'fontweight', 'bold')
    text(800, -550, '(cm yr^{-1})', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold')
    

%% FIGURE 5: BASAL WATER MASKS

    plots                   = {'mask_basal_water_J18_exp' 'mask_basal_water_B19_exp' 'mask_plume_combo_exp' 'mask_basal_water'};
    figure('position', [200 200 1500 800], 'color', 'w', 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; 0.85 0.85 0.85; 1 0.7 0.7; 1 0 0; 0.5 0 0])
    titles                  = {{'(a) Jordan et al. (2018)'; 'basal water identifications'} {'(b) Bowling et al. (2019)'; 'subglacial lakes'} {'(c) Panton and Karlsson (2015) UDRs +'; 'Leysinger-Vieli et al. (2018) basal plumes'} '(d) synthesis'};
    for ii = 1:4
        subplot('position', [(0.005 + ((ii - 1) * 0.25)) 0.11 (0.94 / 4) (0.94 / 1.1)])
        hold on
        plot_tmp            = eval(plots{ii});
        plot_tmp(plot_tmp == 0) ...
                            = 4;
        switch ii
            case {1 4}
                plot_tmp(eval(plots{ii}) >= 1) ...
                            = 5;
                plot_tmp(eval(plots{ii}) >= 5) ...
                            = 6;
                plot_tmp(eval(plots{ii}) >= 10) ...
                            = 7;
            case 2
                plot_tmp(eval(plots{ii}) >= 1) ...
                            = 5;
                plot_tmp(eval(plots{ii}) >= 5) ...
                            = 6;
                plot_tmp(eval(plots{ii}) >= 9) ...
                            = 7;
            case 3
                plot_tmp(eval(plots{ii}) >= 1) ...
                            = 5;
                plot_tmp(eval(plots{ii}) >= 5) ...
                            = 6;
        end
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;        
        imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:num_M19
            plot(M19(jj).X, M19(jj).Y, 'k', 'linewidth', 2)
        end
        plot(x_borehole_frozen, y_borehole_frozen, 'o', 'markersize', 16, 'color', [0 0 0.75], 'linewidth', 2)
        plot(x_borehole_thawed, y_borehole_thawed, 'o', 'markersize', 16, 'color', [0.75 0 0], 'linewidth', 2)
        plot(x_lake, y_lake, 'd', 'markersize', 14, 'color', [0.75 0 0], 'linewidth', 2)        
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        caxis([1 8])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
        if (ii == 3)
            colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', [5 7], 'ytick', 5:7, 'yticklabel', {'' '' ''}, 'ticklength', 0.02, 'fontweight', 'bold')            
        else
            colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', [5 8], 'ytick', 5:8, 'yticklabel', {'' '' '' ''}, 'ticklength', 0.02, 'fontweight', 'bold')
        end
        switch ii
            case 1
                text(1035, -2554, '1-4 / 5-km cell', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
                text(1035, -1655, '5-9 / 5-km cell', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
                text(1035, -777, '10+ / 5-km cell', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
            case 2
                text(1035, -2565, 'low confidence', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
                text(1035, -1550, 'medium confidence', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
                text(1065.5, -1108, {'high or very high'; 'confidence'}, 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270, 'horizontalalignment', 'center')
            case 3
                text(1035, -2231, 'small UDR/plume', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
                text(1035, -915, 'large UDR/plume', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
            case 4
                text(1035, -2555, 'low confidence', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
                text(1035, -1560, 'medium confidence', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
                text(1035, -725, 'high confidence', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
        end
        title(titles{:, ii}, 'fontsize', 22, 'fontweight', 'bold')
    end
    
%% FIGURE 6: SURFACE/DEFORMATION SPEED RATIO

    rb                      = redblue(11, 0.5);
    figure('position', [200 200 1600 800], 'color', 'w', 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; parula(10); rb([1:5 7:11], :)])
    ranges                  = [0 100; 0 100; NaN NaN];
    incs                    = [10 10 NaN];
    plots                   = {'speed_filt' 'speed_def' 'speed_ratio_std_trim'};
    titles                  = {'$\left|\vec{u}_s\right|$' '$u_{def}^T$' '$\left|\vec{u}_s\right| / u_{def}^T$'};
    spaces                  = [9 10 17];
    units                   = {'m yr^{-1}' 'm yr^{-1}' ''};
    [ax, cb]                = deal(NaN(1, 3));
    for ii = 1:3
        ax(ii)              = subplot('position', [(-0.005 + (0.325 * (ii - 1))) 0.02 0.325 0.9]);
        hold on
        if (ii < 3)
            range_tmp       = ranges(ii, 1):incs(ii):ranges(ii, 2);
            plot_tmp        = 3 + discretize(eval(plots{ii}), [-Inf range_tmp(2:(end - 1)) Inf]);
        else
            range_tmp       = [0.5:0.1:1 1.5 2:5];
            plot_tmp        = 13 + discretize(eval(plots{ii}), [-Inf range_tmp(2:(end - 1)) Inf]);
        end
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
        image(x_decim(1, :), y_decim(:, 1), plot_tmp)
        if (ii == 3)
            [~, cu]         = contour(x_decim, y_decim, (speed_filt ./ speed_def), [1 1], 'linewidth', 1.5, 'color', 'w');
        end
        if (ii == 1)
            [~, cuu]        = contour(x_decim, y_decim, (speed_uncert_filt ./ speed_filt), speed_uncert_rel_decay([1 1]), 'linewidth', 1.5, 'color', 'w');
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:num_M19
            plot(M19(jj).X, M19(jj).Y, 'k', 'linewidth', 2)
        end
        plot(x_borehole_frozen, y_borehole_frozen, 'wo', 'markersize', 16, 'markerfacecolor', [0 0 0.75], 'linewidth', 1)
        plot(x_borehole_thawed, y_borehole_thawed, 'wo', 'markersize', 16, 'markerfacecolor', [0.75 0 0], 'linewidth', 1)
        plot(x_lake, y_lake, 'wd', 'markersize', 14, 'markerfacecolor', [0.75 0 0], 'linewidth', 1)
        caxis([1 24])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [])
        if ~isempty(units{ii})
            text(825, -525, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 20)
        end
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-550, -800, repmat(' ', 1, (spaces(ii) + 1)), 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        text(-550, -800, ['(' letters(ii) ')' repmat(' ', 1, spaces(ii))], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold')
        text(-440, -810, titles{ii}, 'color', 'k', 'fontsize', 20, 'interpreter', 'latex')
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 12)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        if (ii == 1)
%             legend(cuu, '$\tilde{u_s} = 25\%$', 'location', 'northeast', 'interpreter', 'latex', 'color', [0.85 0.85 0.85])
        end
        if (ii < 3)
            tick_str        = cell(1, 11);
            for jj = 1:11
                tick_str{jj}= num2str(range_tmp(jj));
            end
            tick_str{end}   = ['>' tick_str{end}];
            cb(ii)          = colorbar('fontsize', 20, 'location', 'eastoutside', 'limits', [4 14], 'ytick', 4:14, 'yticklabel', tick_str, 'ticklength', 0.022);
        else
            cb(ii)          = colorbar('fontsize', 20, 'location', 'eastoutside', 'limits', [14 24], 'ytick', 14:24, 'yticklabel', {'<0.5' '0.6' '0.7' '0.8' '0.9' '1' '1.5' '2' '3' '4' '>5'}, 'ticklength', 0.022);
        end
    end
    annotation('line', [0.9234 0.9335], [0.4688 0.4688], 'color', 'w', 'linewidth', 3)
    pause(1)
    set(cb(1:2), 'limits', [4 14])
    
%% FIGURE 7: BASAL THERMAL STATE AGREEMENT MASKS
    
    tmp1                    = redblue(7);
    tmp                     = [tmp1([2 4:end], :); (2 / 3) 0 0];
    c1                      = [([135 206 235] ./ 255); ([159 89 39] ./ 255); 1 1 1];
    c2                      = [([135 206 235] ./ 255); ([159 89 39] ./ 255); 1 1 1; tmp];
    figure('position', [200 200 1500 800], 'color', 'w', 'renderer', 'zbuffer')
    ax                      = NaN(1, 4);
    titles                  = {'Thawed outlines (standard)' 'Agreement (standard)' 'Agreement (cold bias)' 'Agreement (warm bias)'};
    plots                   = {'' 'mask_agree' 'mask_agree_cold' 'mask_agree_warm'};
    range_tmp               = [-0.25 0 0.25 0.5 0.75 1];
    for ii = 1:4
        ax(ii)              = subplot('position', [(0.02 + ((ii - 1) * 0.245)) 0.11 (0.94 / 4) (0.94 / 1.1)]);
        colormap(c1)
        hold on
        if (ii == 1)
            imagesc(x_decim(1, :), y_decim(:, 1), (mask_combo_decim + 1))
        else
            
            plot_tmp        = 3 + interp1(range_tmp, 1:6, eval(plots{ii}), 'nearest');
            plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
            imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        end
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:num_M19
            plot(M19(jj).X, M19(jj).Y, 'k', 'linewidth', 2)
        end
        if (ii == 1)
            contour(x_decim(1, :), y_decim(:, 1), mask_agree_ismip6, repmat(constraint(3, 2), 1, 2), 'r', 'linewidth', 2);
            contour(x_decim(1, :), y_decim(:, 1), melt_bed_filt, repmat(constraint(1, 2), 1, 2), 'b', 'linewidth', 2);
            contour(x_decim(1, :), y_decim(:, 1), speed_ratio_std, repmat(constraint(2, 2), 1, 2), 'g', 'linewidth', 2);
            contour(x_decim(1, :), y_decim(:, 1), mask_basal_water_std, repmat(constraint(2, 2), 1, 2), 'm', 'linewidth', 2);
            psr             = plot(NaN, NaN, 'r', 'linewidth', 2);
            psf             = plot(NaN, NaN, 'b', 'linewidth', 2);
            pud             = plot(NaN, NaN, 'g', 'linewidth', 2);
            pmg             = plot(NaN, NaN, 'm', 'linewidth', 2);
        end
        pbf                 = plot(x_borehole_frozen, y_borehole_frozen, 'ko', 'markersize', 16, 'markerfacecolor', [0 0 0.75], 'linewidth', 1);
        pbt                 = plot(x_borehole_thawed, y_borehole_thawed, 'ko', 'markersize', 16, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
        psl                 = plot(x_lake, y_lake, 'kd', 'markersize', 14, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
        if (ii > 1)
            set([pbf pbt psl], 'color', 'w')
        end
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
        title(['(' letters(ii) ') ' titles{ii}], 'color', 'k', 'fontweight', 'bold')
        switch ii
            case 1
                caxis([1 3])
            otherwise
                caxis([1 10])
                colormap(ax(ii), c2)
        end
        if (ii == 2)
            lg              = legend([pbf pbt psl psr psf pud pmg], 'borehole (frozen bed)', 'borehole (thawed bed)', 'subglacial lake', 'ISMIP6 thawed', '$\dot{m}$ = 1 cm yr$^{-1}$', '$\gamma_{min} = 1$', 'basal water', 'interpreter', 'latex');
            set(lg, 'position', [0.0676 0.0142 0.1454 0.1991])
            cb              = colorbar('fontsize', 20, 'location', 'southoutside', 'xlim', [4 10], 'xtick', 4:10, 'xticklabel', {},'position', [0.265 0.08 (0.76 + (0.91 / 4) - 0.26) 0.025], 'ticklength', 0.018);
            tick_str        = {'-1 frozen' '0 no agreement' '+1 thawed' '+2 thawed' ' +3 thawed' '+4 all agree thawed'};
            x_tick          = [-450 200 1100 1800 2550 3150];
            for jj = 1:6
                text(x_tick(jj), -3600, tick_str{jj}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold')
            end
        end
    end
    linkaxes(ax)
    pause(1)
    set(ax(2:3), 'colormap', get(ax(4), 'colormap'))
    
%% FIGURE 8: GBaTSv2 LIKELY BASAL THERMAL STATE
    
    rb                      = redblue(3);
    rb(2, :)                = [0.8 0.8 0.8];
    rb2                     = [0.5 0.5 1; 0.9 0.9 0.9; 1 0.5 0.5];
    plots                   = {'GBaTSv1.mask_likely_filled' 'mask_likely_filled' '(mask_likely_filled - GBaTSv1.mask_likely_filled)'};
    titles                  = {'GBaTSv1 (M16)' 'GBaTSv2 (this study)' 'difference (GBaTSv2 - GBaTSv1)'};
    figure('position', [200 200 1600 800], 'color', 'w', 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; rb; rb2])
    for ii = 1:3
        subplot('position', [(-0.005 + (0.325 * (ii - 1))) 0.02 0.325 0.9])
        hold on
        if (ii < 3)
            plot_tmp        = interp1(-1:1, 4:6, eval(plots{ii}), 'nearest', 'extrap');
        else
            plot_tmp        = interp1(-1:1, 7:9, eval(plots{ii}), 'nearest', 'extrap');
        end
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
        imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:num_M19
            plot(M19(jj).X, M19(jj).Y, 'k', 'linewidth', 2)
        end
        if (ii == 1)
            plot([x_borehole_frozen(1:8) x_borehole_thawed(10)], [y_borehole_frozen(1:8) y_borehole_thawed(10)], 'wo', 'markersize', 16, 'markerfacecolor', [0 0 0.75], 'linewidth', 1)
            plot(x_borehole_thawed(1:8), y_borehole_thawed(1:8), 'wo', 'markersize', 16, 'markerfacecolor', [0.75 0 0], 'linewidth', 1)
            plot(x_lake(1:4), y_lake(1:4), 'wd', 'markersize', 14, 'markerfacecolor', [0.75 0 0], 'linewidth', 1)
        else
            plot(x_borehole_frozen, y_borehole_frozen, 'wo', 'markersize', 16, 'markerfacecolor', [0 0 0.75], 'linewidth', 1)
            plot(x_borehole_thawed, y_borehole_thawed, 'wo', 'markersize', 16, 'markerfacecolor', [0.75 0 0], 'linewidth', 1)
            plot(x_lake, y_lake, 'wd', 'markersize', 14, 'markerfacecolor', [0.75 0 0], 'linewidth', 1)
        end
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 13)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        caxis([1 10])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
        if (ii < 3)
            colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', [4 7], 'ytick', 4:7, 'yticklabel', {'' '' '' ''}, 'ticklength', 0.02, 'fontweight', 'bold')
            text(1050, -900, 'likely thawed', 'fontsize', 22, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
            text(1050, -1850, 'uncertain', 'fontsize', 22, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
            text(1050, -2700, 'likely frozen', 'fontsize', 22, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
        else
            colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', [7 10], 'ytick', 7:10, 'yticklabel', {'' '' '' ''}, 'ticklength', 0.02, 'fontweight', 'bold')
            text(1050, -750, 'more likely thawed', 'fontsize', 22, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
            text(1050, -1800, 'no change', 'fontsize', 22, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
            text(1050, -2550, 'more likely frozen', 'fontsize', 22, 'fontweight', 'bold', 'color', 'k', 'rotation', 270)
        end
        if (ii == 1)
            text(192, -1945, 'S', 'color', 'w', 'fontsize', 18, 'fontweight', 'bold')
            text(-27, -1683, 'NG', 'color', 'w', 'fontsize', 18, 'fontweight', 'bold')
        end
        title(['(' letters(ii) ') ' titles{ii}], 'fontweight', 'bold')
    end

%% FIGURE 9: ISMIP6 VS SEARISE AGREEMENT
    
    tmp                     = redblue(12, 0.5);
    tmp2                    = flipud(hot(10));
    figure('position', [200 200 1600 800], 'color', 'w', 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; tmp([1:5 8:12], :); tmp2])
    titles                  = {{'(a) Change in ice thickness'; 'from SeaRISE to ISMIP6 (m)'} {'(b) Change in agreement that bed is likely thawed'; 'from SeaRISE to ISMIP6 (%)'}};
    for ii = 1:2
        subplot('position', [(-0.005 + (0.35 * (ii - 1))) 0.02 0.325 0.9])
        hold on
        switch ii
            case 1
                range_tmp   = linspace(-500, 500, 11);
                plot_tmp    = NaN(num_decim_y, num_decim_x);
                plot_tmp(~isnan(SeaRISE.thick)) ...
                            = discretize((BM3.thick_interp(~isnan(SeaRISE.thick)) - SeaRISE.thick(~isnan(SeaRISE.thick))), [-Inf range_tmp(2:(end - 1)) Inf]);
            case 2
                range_tmp   = linspace(-1, 1, 11);
                plot_tmp    = NaN(num_decim_y, num_decim_x);
                plot_tmp(mask_agree_diff < 0) ...
                            = discretize(mask_agree_diff(mask_agree_diff < 0), range_tmp, 'includededge', 'left');
                plot_tmp(mask_agree_diff == 0) ...
                            = find(range_tmp == 0);
                plot_tmp(mask_agree_diff > 0) ...
                            = discretize(mask_agree_diff(mask_agree_diff > 0), range_tmp, 'includededge', 'right');
        end
        plot_tmp            = 3 + plot_tmp;
        plot_tmp(isnan(plot_tmp) | ~mask_gris_decim) ...
                            = mask_combo_decim(isnan(plot_tmp) | ~mask_gris_decim) + 1;
        imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:num_M19
            plot(M19(jj).X, M19(jj).Y, 'k', 'linewidth', 2)
        end
        pbf                 = plot(x_borehole_frozen, y_borehole_frozen, 'wo', 'markersize', 16, 'markerfacecolor', [0 0 0.75], 'linewidth', 1);
        pbt                 = plot(x_borehole_thawed, y_borehole_thawed, 'wo', 'markersize', 16, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
        psl                 = plot(x_lake, y_lake, 'wd', 'markersize', 14, 'markerfacecolor', [0.75 0 0], 'linewidth', 1);
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 13)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        caxis([1 24])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
        tick_str            = cell(1, 11);
        switch ii
            case 1
                for jj = 1:11
                    tick_str{jj}= num2str(range_tmp(jj));
                end
            case 2
                for jj = 1:11
                    tick_str{jj}= num2str(1e2 * range_tmp(jj));
                end
        end
        tick_str{1}         = ['<' tick_str{1}];
        tick_str{end}       = ['>' tick_str{end}];
        colorbar('fontsize', 20, 'location', 'eastoutside', 'limits', [4 14], 'ytick', 4:14, 'ticklength', 0.021, 'fontweight', 'bold', 'yticklabel', tick_str)
        title(titles{ii}, 'fontsize', 22, 'fontweight', 'bold', 'color', 'k', 'horizontalalignment', 'center')
    end
    subplot('position', [0.70 0.2 0.30 0.6])
    hold on
    range_tmp               = linspace(0, 0.5, 11);
    p                       = imagesc('xdata', -2000:100:2000, 'ydata', -150:10:150, 'cdata', (14 + discretize(model_change_bin, [-Inf range_tmp(2:(end - 1)) Inf])));
    set(p, 'alphadata', (model_change_bin > 0))
    axis square
    axis([-2000 2000 -150 150])
    caxis([1 24])
    tick_str                = cell(1, 11);
    for ii = 1:2:11
        tick_str{ii}        = num2str(range_tmp(ii));
    end
    tick_str{end}           = ['>' tick_str{end}];
    colorbar('fontsize', 20, 'location', 'southoutside', 'limits', [15 25], 'ytick', 15:25, 'ticklength', 0.04, 'fontweight', 'bold', 'yticklabel', tick_str)
    set(gca, 'fontsize', 20, 'fontweight', 'bold', 'layer', 'top')
    xlabel('Change in thickness (m)')
    ylabel('Change in agreement on thawed bed (%)')
    title('(c)', 'fontsize', 22, 'fontweight', 'bold')
    text(-1400, -265, 'Fraction of grid cells (%)', 'fontsize', 22, 'fontweight', 'bold')
    box on
    grid on
        
end