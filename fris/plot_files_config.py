#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# Centralized line and time property definitions
Lwides1 = 2.0
Lwides6 = 1.5

lstyl1 = '--'
lstyl6 = '-'

S6_start_yr = 0

files_config = [
    {
        'filename': 'bulk_tseries_2D_1by1_F8_Spin1.nc',
        'label': 'F8 (Spin1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'orange',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False,
        'meshres': 8.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F8_Spin6.nc',
        'label': 'F8 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'brown',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True,
        'meshres': 8.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F4_Spin1.nc',
        'label': 'F4 (Spin1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'lightskyblue',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False,
        'meshres': 4.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F4_Spin6.nc',
        'label': 'F4 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'royalblue',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True,
        'meshres': 4.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F2_Spin6.nc',
        'label': 'F2 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'forestgreen',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True,
        'meshres': 2.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F2_Spin1p1.nc',
        'label': 'F2 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'yellowgreen',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False,
        'meshres': 2.0
    },
{
        'filename': 'bulk_tseries_2D_1by1_F2_Spin1p2.nc',
        'label': 'F2 (Spin1p2)',
        'start_year': 2,
        'start_month': 11,
        'color': 'yellowgreen',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False,
        'meshres': 2.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F1_Spin1p1.nc',
        'label': 'F1 (Spin1p1)',
        'start_year': 0,
        'start_month': 1,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False,
        'meshres': 1.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F1_Spin1p2.nc',
        'label': 'F1 (Spin1p2)',
        'start_year': 0,
        'start_month': 12,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False,
        'meshres': 1.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F1_Spin1p3.nc',
        'label': 'F1 (Spin1p3)',
        'start_year': 4,
        'start_month': 6,
        'color': 'darkgray',
        'linewidth': Lwides1,
        'linestyle': lstyl1,
        'in_legend': False,
        'meshres': 1.0
    },
    {
        'filename': 'bulk_tseries_2D_1by1_F1_Spin6.nc',
        'label': 'F1 (Spin6)',
        'start_year': S6_start_yr,
        'start_month': 1,
        'color': 'black',
        'linewidth': Lwides6,
        'linestyle': lstyl6,
        'in_legend': True,
        'meshres': 1.0
    },
]



if __name__ == "__main__":

    for config in files_config:
        print(config['filename'])




