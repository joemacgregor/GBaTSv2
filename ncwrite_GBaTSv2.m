% NCWRITE_GBaTSv2 Make NetCDF database of Greenland basal thermal state v2.
% 
% Joe MacGregor
% Last updated: 18 October 2021

clear

path_nc                     = '/Users/jamacgre/Documents/research/manuscripts/GBaTSv2/misc/';
file_nc                     = 'GBaTSv2.nc';
nc_ver                      = 0.3;
val_nodata                  = NaN;

load ~/Documents/data/mat/bed2 mask_agree mask_agree_cold mask_agree_warm mask_agree_ismip6 mask_agree_ismip6_cold mask_agree_ismip6_warm mask_basal_water mask_likely_filled mask_D1_decim melt_bed_filt melt_bed_lo_filt melt_bed_hi_filt ...
                               num_decim_x num_decim_y speed_ratio_std_trim speed_ratio_min_trim speed_ratio_max_trim x_decim y_decim

mapping_bm                  = ncinfo('/Users/jamacgre/Documents/research/data/greenland/BedMachine/BedMachineGreenland-2021-04-20.nc');
mapping_bm                  = mapping_bm.Variables(1);

[x_decim, y_decim]          = deal((1e3 .* x_decim(1, :))', (flipud((1e3 .* y_decim(:, 1))))); % km to m

% ISMIP6 agreement mask (Figure 3)
mask_agree_ismip6           = rot90((10 .* mask_agree_ismip6), -1);
mask_agree_ismip6_cold      = rot90((10 .* mask_agree_ismip6_cold), -1);
mask_agree_ismip6_warm      = rot90((10 .* mask_agree_ismip6_warm), -1);

% basal melt rate from radiostratigraphy (Figure 4)
[melt_bed_filt(~mask_D1_decim), melt_bed_lo_filt(~mask_D1_decim), melt_bed_hi_filt(~mask_D1_decim)] ...
                            = deal(NaN);
melt_bed                    = rot90(melt_bed_filt, -1);
melt_bed_min                = rot90(melt_bed_lo_filt, -1);
melt_bed_max                = rot90(melt_bed_hi_filt, -1);
clear melt_bed_*_filt

% basal water mask (Figure 5d)
mask_basal_water            = rot90(mask_basal_water, -1);
clear mask_basal_water_*

% ratio of surface speed to SIA temperate column deformation speed (Figure 6c)
speed_ratio_std_trim        = rot90(speed_ratio_std_trim, -1);
speed_ratio_min_trim        = rot90(speed_ratio_min_trim, -1);
speed_ratio_max_trim        = rot90(speed_ratio_max_trim, -1);

% agreement mask (Figure 7)
mask_agree                  = rot90((4 .* mask_agree), -1);
mask_agree_cold             = rot90((4 .* mask_agree_cold), -1);
mask_agree_warm             = rot90((4 .* mask_agree_warm), -1);

% likely basal thermal state (Figure 8)
mask_likely_filled          = rot90(mask_likely_filled, -1);

%%

if exist([path_nc file_nc], 'file')
    eval(['!mv ' path_nc file_nc ' ' path_nc file_nc(1:(end - 3)) '_old.nc'])
end

% variable writing to NetCDF file
nccreate([path_nc file_nc], 'mapping', 'Datatype', 'char', 'Format', 'netcdf4')
nccreate([path_nc file_nc], 'x', 'Dimensions', {'x' num_decim_x}, 'Datatype', 'double')
ncwrite([path_nc file_nc], 'x', x_decim)
nccreate([path_nc file_nc], 'y', 'Dimensions', {'y' num_decim_y}, 'Datatype', 'double')
ncwrite([path_nc file_nc], 'y', y_decim)

nc_vars                     = {'agreement_ismip6'                   'mask_agree_ismip6'         'ISMIP6 standard agreement mask (threshold: -1.0 degrees Celsius); sum of model agreement, where for each model a thawed basal thermal state is +1 and negative is -1';
                               'agreement_ismip6_cold'              'mask_agree_ismip6_cold'    'ISMIP6 cold-bias agreement mask (threshold: -0.5 degrees Celsius); sum of model agreement, where for each model a thawed basal thermal state is +1 and negative is -1';
                               'agreement_ismip6_warm'              'mask_agree_ismip6_warm'    'ISMIP6 warm-bias agreement mask (threshold: -1.5 degrees Celsius); sum of model agreement, where for each model a thawed basal thermal state is +1 and negative is -1';
                               'basal_melt'                         'melt_bed'                  'standard basal melt rate inferred from Nye+melt 1-D model';
                               'basal_melt_min'                     'melt_bed_min'              'minimum basal melt rate inferred from Nye+melt 1-D model';
                               'basal_melt_max'                     'melt_bed_max'              'maximum basal melt rate inferred from Nye+melt 1-D model';
                               'basal_water'                        'mask_basal_water'          'basal water mask that synthesizes several published inferences of basal water; count of basal water identifications';
                               'speed_ratio'                        'speed_ratio_std_trim'      'standard ratio of observed surface speed to modeled deformation of a temperate ice column; dimensionless ratio';
                               'speed_ratio_min'                    'speed_ratio_min_trim'      'minimum ratio of observed surface speed to modeled deformation of a temperate ice column; dimensionless ratio';
                               'speed_ratio_max'                    'speed_ratio_max_trim'      'maximum ratio of observed surface speed to modeled deformation of a temperate ice column; dimensionless ratio';
                               'agreement_basal_thermal_state'      'mask_agree'                'standard four-method agreement mask; total balance of methods indicating a specific basal thermal state (negative:frozen, positive:thawed';
                               'agreement_basal_thermal_state_cold' 'mask_agree_cold'           'cold-bias four-method agreement mask; total balance of methods indicating a specific basal thermal state (negative:frozen, positive:thawed';
                               'agreement_basal_thermal_state_warm' 'mask_agree_warm'           'warm-bias four-method agreement mask; total balance of methods indicating a specific basal thermal state (negative:frozen, positive:thawed';
                               'likely_basal_thermal_state'         'mask_likely_filled'        'likely basal thermal state; -1: likely frozen; 0: uncertain; +1: likely thawed'};
for ii = 1:size(nc_vars, 1)
    nccreate([path_nc file_nc], nc_vars{ii, 1}, 'Dimensions', {'x' num_decim_x 'y' num_decim_y}, 'Datatype', 'double', 'FillValue', val_nodata)
    ncwrite([path_nc file_nc], nc_vars{ii, 1}, eval(nc_vars{ii, 2}))
end

% global attributes
ncwriteatt([path_nc file_nc], '/', 'Conventions', 'CF-1.7')
ncwriteatt([path_nc file_nc], '/', 'Title', 'GBaTSv2: A revised synthesis of the likely basal thermal state of the Greenland Ice Sheet')
ncwriteatt([path_nc file_nc], '/', 'Author', 'Joseph MacGregor')
ncwriteatt([path_nc file_nc], '/', 'Version', [datestr(now) ' (v' num2str(nc_ver) ')'])
ncwriteatt([path_nc file_nc], '/', 'Data_citation', 'MacGregor et al., in prep.')
ncwriteatt([path_nc file_nc], '/', 'nx', num_decim_x)
ncwriteatt([path_nc file_nc], '/', 'ny', num_decim_y)
ncwriteatt([path_nc file_nc], '/', 'Projection', 'EPSG:3413 NSIDC Sea Ice Polar Stereographic North, centered on 70N, 45W')
ncwriteatt([path_nc file_nc], '/', 'proj4', '+init=epsg:3413')
ncwriteatt([path_nc file_nc], '/', 'xmin', min(x_decim))
ncwriteatt([path_nc file_nc], '/', 'ymax', max(y_decim))
ncwriteatt([path_nc file_nc], '/', 'spacing', 5e3)
ncwriteatt([path_nc file_nc], '/', 'no_data', val_nodata)

% variable-specific attributes
for ii = 2:length(mapping_bm.Attributes) % skip geoid
    ncwriteatt([path_nc file_nc], 'mapping', mapping_bm.Attributes(ii).Name, mapping_bm.Attributes(ii).Value)
end
ncwriteatt([path_nc file_nc], 'mapping', 'long_name', 'mapping')
ncwriteatt([path_nc file_nc], 'mapping', 'spatial_ref', 'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",SOUTH],AXIS["Northing",SOUTH],AUTHORITY["EPSG","3413"]]')
ncwriteatt([path_nc file_nc], 'mapping', 'crs_wkt', 'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",SOUTH],AXIS["Northing",SOUTH],AUTHORITY["EPSG","3413"]]')
ncwriteatt([path_nc file_nc], 'mapping', 'proj4text', '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
ncwriteatt([path_nc file_nc], 'mapping', 'GeoTransform', [-632500 5000 0 -667500 0 -5000])
ncwriteatt([path_nc file_nc], 'x', 'long_name', 'Cartesian x-coordinate')
ncwriteatt([path_nc file_nc], 'x', 'standard_name', 'projection_x_coordinate')
ncwriteatt([path_nc file_nc], 'x', 'units', 'meter')
ncwriteatt([path_nc file_nc], 'y', 'long_name', 'Cartesian y-coordinate')
ncwriteatt([path_nc file_nc], 'y', 'standard_name', 'projection_y_coordinate')
ncwriteatt([path_nc file_nc], 'y', 'units', 'meter')
for ii = 1:size(nc_vars, 1)
    ncwriteatt([path_nc file_nc], nc_vars{ii, 1}, 'grid_mapping', 'mapping')
    ncwriteatt([path_nc file_nc], nc_vars{ii, 1}, 'long_name', nc_vars{ii, 1})
    if contains(nc_vars{ii, 1}, 'basal_melt')
        ncwriteatt([path_nc file_nc], nc_vars{ii, 1}, 'units', 'meters per year')
    end
    ncwriteatt([path_nc file_nc], nc_vars{ii, 1}, 'description', nc_vars{ii, 3})
end