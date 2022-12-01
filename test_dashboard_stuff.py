import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import pandas as pd
import plotly.express as px
from jupyter_dash import JupyterDash

import plotly.graph_objects as go

import csv
import h5py
import geopandas as gpd
import astropy.time as astro_time

import time

from wmm import WMM

countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options



class OAT:

    def __init__(self, fn, path = 'save_sim/', output_path = 'outputs/'):

        print('Initializing OAT simulation..')

        # set the font globally
        #plt.rcParams.update({'font.family':'sans-serif'})

        self.sim_path = path
        self.out_path = output_path

        f = h5py.File(self.sim_path + fn, 'r')
        print('Loading ' + fn + '..')

        self.dt = f.attrs['dt']

        self.X = f['states']['ECI']['X']
        self.Y = f['states']['ECI']['Y']
        self.Z = f['states']['ECI']['Z']
        self.LALN = f['states']['angular']['LALN']
        self.H = f['states']['ECI']['H']
        self.OE_ = f['states']['angular']['OE']

        self.B = f['B']['B']
        self.Bx = f['B']['Bx']
        self.By = f['B']['By']
        self.Bz = f['B']['Bz']

        self.times_jd = astro_time.Time(f['times']['JD'], format = 'jd').jd

        self.times_utc = astro_time.Time(f['times']['JD'], format = 'jd').datetime

        print('data successfully loaded..')

oat = OAT('test.hdf5')

print(oat.LALN)

print(oat.LALN[0])
print(oat.LALN[0][1])
print(oat.LALN[0,1])

'''for something in oat.LALN:
    print(something)'''