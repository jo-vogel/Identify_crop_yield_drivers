#!/usr/bin/env python

"""
This script reads daily data from three gridpoints and their annual crop
yield and produces Figures 2, B2, B3 in the paper.

INPUT:
- meteo_daily_FR_ensemblemode.nc
- meteo_daily_CH2_ensemblemode.nc
- meteo_daily_US3_ensemblemode.nc
- crop_yield_NH.nc

OUTPUT:
Figures 2, B2, B3 in the paper.

author: Christoph A Sauter
email: christoph.sauter@strath.ac.uk
date: 22. Sep 2020
Python version: 3.7
"""

from netCDF4 import Dataset
import numpy as np
from numba import jit
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
    


# -------------------------
# For reading the input
def read_nc_data(data_dir, file_name, var_list):
    """
    reads netCDF data from a netCDF file into arrays
    :param data_dir: directory of the file
    :param file_name: name of netCDF file
    :param var_list: name of variables in a list, e.g. ['temp', 'precip',]
    :return: longitude, latitude, time, variables as entries in a dict
    """
    ncdata = Dataset(data_dir + file_name, mode='r')
    var_dict = {}
    for var in var_list:
        var_dict[var] = ncdata.variables[var][:]
    ncdata.close()
    return var_dict


# -------------------------
# Functions for general calculations
def get_extremes_1gp(met_var, crop, threshold):
    # divide into low yield and rest
    neg_extreme_indices = crop[:] < np.quantile(crop[:], threshold)
    avr_and_pos_indices = ~np.array(neg_extreme_indices)
    met_var_neg = met_var[neg_extreme_indices]
    met_var_av = met_var[avr_and_pos_indices]
    return met_var_neg, met_var_av


@jit(nopython=True)
def running_sum(x, running_days):
    # Creates a running sum of a timeseries.
    rs = np.zeros(x.shape)
    rs[:, :] = np.nan
    for year in range(0, x.shape[0]):
        for day in range(running_days, x.shape[1]):
            rs[year, day] = np.sum(x[year, day-running_days:day])
    return rs


# -------------------------
# Functions for plotting
def subplot_country_location(ax0, loc_lon, loc_lat):
    # Creates a map with dot for gridpoint
    ax0.set_extent([loc_lon - 32, loc_lon + 32, loc_lat - 20, loc_lat + 20],
                   crs=ccrs.PlateCarree())
    ax0.coastlines(resolution='50m')
    plt.scatter(loc_lon, loc_lat, c='orangered', s=70)
    ax0.set_xticks(np.arange(np.ceil((loc_lon - 32) / 10) * 10,
                             np.ceil((loc_lon + 32) / 10) * 10, 15),
                   crs=ccrs.PlateCarree())
    ax0.set_yticks(np.arange(np.ceil((loc_lat - 20) / 10) * 10,
                             np.ceil((loc_lat + 20) / 10) * 10, 10),
                   crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(degree_symbol=u'\u00B0 ')
    lat_formatter = LatitudeFormatter(degree_symbol=u'\u00B0 ')
    ax0.xaxis.set_major_formatter(lon_formatter)
    ax0.yaxis.set_major_formatter(lat_formatter)
    ax0.xlabel_style = {'size': 1}
    return None


def plot_mean_with_percentiles(x, ax, alpha, plot_color, percentile_value, 
                               start):
    """
    :param x: variable with axis 0:
    :param alpha: alpha value of percentile shading
    :param plot_color: plot color
    :param percentile_value: draws shades from upper to lower percentile
    :param start: sowing date day (days since 01.01.2035)
    :return: the figure
    """

    x_lower_bound = np.percentile(x, percentile_value, axis=0)
    x_upper_bound = np.percentile(x, 100-percentile_value, axis=0)

    months = mdates.MonthLocator()
    months_fmt = mdates.DateFormatter('%b')
    # Function for displaying only first letter of month
    def m_fmt(x, pos=None):
        return months_fmt(x)[0]
    
    base_date = datetime.strptime('2035-01-01', '%Y-%m-%d')
    start_date = base_date + timedelta(days=start)
    end_date = start_date + timedelta(days=x.shape[1])
    date_range = pd.date_range(start_date, end_date)
    
    
    ax.plot(date_range[:-1], np.mean(x, axis=0), 
             linewidth=2, color=plot_color)

    ax.fill_between(date_range[:-1], x_lower_bound, x_upper_bound,
                     alpha=alpha, color=plot_color, linewidth=0,
                     label='_nolegend_')

    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_major_formatter(FuncFormatter(m_fmt))
    ax.set_xlim(start_date - timedelta(days=10), 
                end_date + timedelta(days=10))
    
    for pos in ['right', 'top']:
        plt.gca().spines[pos].set_visible(False)
    
    return None



# -------------------------
# BEGIN OF SCRIPT
# -------------------------

# EDIT TO DIRECTORY WHERE THE FILES ARE:
data_dir = ('/Users/christoph/Desktop/PhD/DAMOCLES_training_school/'
            'WorkingGroup1/nc_data/')

met_var_list = ['dps', 'pr', 'rsds', 'sfcwind', 'tasmax', 'tasmin', 'vpd']

crop = read_nc_data(data_dir, file_name='crop_yield_NH.nc',
                    var_list=['yield', 'lat', 'lon'])

gp = ['FR', 'CH2', 'US3']
# France 47.7 ​N, 1.1 ​E
# China 30.8 ​N, 118.1 E
# USA  44.3 ​N, 90.0 W

data_gp = {}
for i, loc in enumerate(gp):
    fname = 'meteo_daily_' + gp[i] + '_ensemblemode.nc'
    data_gp[loc] = read_nc_data(data_dir, file_name=fname,
                                var_list=met_var_list)

gp_locs_fullnames = {'FR': 'France',
                     'CH2': 'China',
                     'US3': 'USA'}

# Grid point locations (as index)
gp_locs = {'FR': [42, 160],
           'CH2': [27, 264],
           'US3': [39, 79]}

# Sowing date; Beginning of growing season
gp_start = {'FR': 304,
            'CH2': 272,
            'US3': 263}

# Length of the longest growing season
gp_length_max = {'FR': 310,
                 'CH2': 258,
                 'US3': 333}


# divide all data into negative extremes and others (threshold 5th percentile)
th_percentile = 0.05
data_extr = {}
for loc in gp_locs:
    first_loop = True
    for met_var in met_var_list:
        var_neg, var_av = get_extremes_1gp(data_gp[loc][met_var],
                                           crop['yield'][:, gp_locs[loc][0],
                                           gp_locs[loc][1]],
                                           threshold=th_percentile)
        if first_loop is True:
            data_extr[loc] = {met_var+'_neg': var_neg,
                              met_var+'_av': var_av}
            first_loop = False
        else:
            data_extr[loc].update({met_var+'_neg': var_neg,
                                   met_var+'_av': var_av})

    # Convert from Kelvin to Celsius
    temp_vars = ['dps_neg', 'dps_av', 
                 'tasmax_neg', 'tasmax_av', 
                 'tasmin_neg', 'tasmin_av']
    
    for tvar in temp_vars:
        data_extr[loc][tvar] -= 273.16


# PLOT SETTINGS
orange = '#d95f02'
blue = '#7570b3'
alpha = 0.3
var_names_plot = {'dps': 'T$_d$  [$^\circ$C]',
                  'pr_rsum': 'Pr, Running sum: 30 days  [mm]',
                  'rsds': 'Rad  [W m$^{-2}$]',
                  'sfcwind': 'Wind  [m/s]',
                  'tasmax': r'T$_{max}$  [$^\circ$C]',
                  'tasmin': r'T$_{min}$  [$^\circ$C]',
                  'vpd': 'VPD  [hPa]'}
alphabet = 'abcdefghijklmnopqrstuvwxyz'


# PLOT YEARLY COMPOSITES
# Plot yearly evolution of each LY gridpoint for all variables
# November: day [:, 304:304+365]
lab_fontsize = 22
ax_fontsize = 18
var_lost_plot = ['dps', 'pr_rsum', 'rsds', 'sfcwind', 'tasmax', 'tasmin',
                 'vpd']

for country in gp:
    print('Plotting grid point in ' + gp_locs_fullnames[country])

    fig = plt.figure(figsize=(12*1.1, 16*1.1))

    var_counter = 0
    data_extr[country]['pr_rsum_neg'] = running_sum(data_extr[
                                                    country]['pr_neg'],
                                                    running_days=30)
    data_extr[country]['pr_rsum_av'] = running_sum(data_extr[
                                                   country]['pr_av'],
                                                   running_days=30)
    loc_lon = crop['lon'][gp_locs[country][1]]
    loc_lat = crop['lat'][gp_locs[country][0]]

    for ax_v in range(4):
        for ax_h in range(2):
            if var_counter == 0:
                # Plot map with gridpoint location
                ax0 = fig.add_subplot('42'+str(var_counter+1),
                                      projection=ccrs.PlateCarree())
                ax0.tick_params(labelsize=ax_fontsize)
                plt.title('(' + alphabet[var_counter] + ')',
                          loc='left', pad=10, fontsize=lab_fontsize)
                subplot_country_location(ax0, loc_lon, loc_lat)
            else:
                # Plot variables
                var = var_lost_plot[var_counter - 1]
                ax = fig.add_subplot('42' + str(var_counter+1))
                ax.tick_params(labelsize=ax_fontsize)
                plot_mean_with_percentiles(data_extr[country][var
                    + '_neg'][:, gp_start[country]:gp_start[country]
                                 +gp_length_max[country]], 
                    ax=ax, alpha=alpha, plot_color=orange,
                    percentile_value=10, start=gp_start[country],
                    )
                plot_mean_with_percentiles(data_extr[country][var
                    + '_av'][:, gp_start[country]:gp_start[country]
                                +gp_length_max[country]], 
                    ax=ax, alpha=alpha, plot_color=blue, percentile_value=10, 
                    start=gp_start[country]
                    )
                plt.ylabel(var_names_plot[var], fontsize=ax_fontsize+2)
                plt.title('(' + alphabet[var_counter] + ')',
                          loc='left', fontsize=lab_fontsize)
                if var_counter == 1:
                    plt.legend(('Bad years', 'Normal years'),
                               frameon=False, loc=2, fontsize=ax_fontsize)
            var_counter = var_counter + 1
    fig.canvas.draw()  # Work-around for when tight_layout() produces an error
    plt.tight_layout()
    # fig.savefig(data_dir[:-8] + 'Plots/' + country +
    #             '_yearly_composite_rev.pdf')
    plt.show()
       
