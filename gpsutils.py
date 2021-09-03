import os, csv
import re
import matplotlib.pyplot as plt
import numpy as np
from dateutil.parser import isoparse
from utils import return_all_filenames_in_path
from Strain_2D.Strain_Tools.strain.utilities import GPSData
from numpy.linalg import inv

class GPSNetwork(object):

    def __init__(self, filenames, prefix = '.'):
        self._network_data = {}
        for filename in filenames:
            name, lat, lon = None, None, None
            try:
                with open(os.path.join(prefix, filename), 'r') as file:
                    reader = csv.reader(file, delimiter = ',')
                    for row in reader:
                        line = row[0]
                        match = re.search('Latitude: ([+-]?([0-9]*[.])?[0-9]+)', line)
                        if match:
                            lat = float(match.group(1))
                        match = re.search('Longitude: ([+-]?([0-9]*[.])?[0-9]+)', line)
                        if match:
                            lon = float(match.group(1))
                        match = re.search('Source File: ([A-Z0-9]{4})', line)
                        if match:
                            name = match.group(1)
                        if lat is not None and lon is not None and name is not None:
                            break
            except:
                print('Problem with ', filename)
                break
            if name is not None:
                self._network_data[name] = [lon, lat]
            else:
                print('Problem with {name}'.format(name = filename))

    def plot(self, axis = None, color = 'k.', textcolor = 'k'):
        for key, value in self._network_data.items():
            if axis is None or (value[0] > axis[0] and value[0] < axis[1] and value[1] > axis[2] and value[1] < axis[3]):
                plt.plot(value[0], value[1], color)
                plt.text(value[0], value[1], key)
                plt.axis('equal')


    def get_lat_lon_for_name(self, name):
        return self._network_data.get(name, [np.nan, np.nan])

    def get_stations(self):
        from Strain_2D.Strain_Tools.strain.utilities import Stations
        return [Stations(elon = value[0], nlat = value[1], name = key) for key,value in self._network_data.items()]

    def get_strain_calculator(self, grid_inc, strain_range = None, strain_method = 'delaunay', **kwargs):
        from Strain_2D.Strain_Tools.strain.configure_functions import Params
        stations = self.get_stations()
        if strain_range is None:
            data_range = self.data_range
            xrange = (data_range[1]-data_range[0])*1.1
            yrange = (data_range[3]-data_range[2])*1.1
            meanx = (data_range[1] + data_range[0]) / 2.0
            meany = (data_range[3] + data_range[2]) / 2.0
            strain_range = [meanx-xrange/2.0, meanx+xrange/2.0,meany-yrange/2.0,meany+yrange/2.0]
        myParams = Params(strain_method=strain_method, input_file=None, range_strain= strain_range, range_data = self.data_range, inc = grid_inc, outdir = None, method_specific=kwargs)
        strain_calculator = None
        if strain_method == 'delaunay_flat':
            from Strain_2D.Strain_Tools.strain.models.strain_delaunay_flat import delaunay_flat
            strain_calculator = delaunay_flat(myParams)
        elif strain_method == 'delaunay':
            from Strain_2D.Strain_Tools.strain.models.strain_delaunay import delaunay
            strain_calculator = delaunay(myParams)
        elif strain_method == 'geostats':
            from Strain_2D.Strain_Tools.strain.models.strain_geostats import geostats
            strain_calculator = geostats(myParams)
        elif strain_method == 'huang':
            from Strain_2D.Strain_Tools.strain.models.strain_huang import huang
            strain_calculator = huang(myParams)
        elif strain_method == 'gpsgridder':
            from Strain_2D.Strain_Tools.strain.models.strain_gpsgridder import gpsgridder
            strain_calculator = gpsgridder(myParams)
        else:
            raise Exception('Invalid strain method.')
        strain_calculator.configure_network(stations)
        return strain_calculator

    @property
    def data_range(self):
        lons = np.array([item[0] for key,item in self._network_data.items()])
        lats = np.array([item[1] for key,item in self._network_data.items()])
        return [np.min(lons), np.max(lons), np.min(lats), np.max(lats)]

    @property
    def yrange(self):
        return self._yrange

class GPSTimeSeries(object):

    def __init__(self, filename, prefix = '.'):
        self._name = None
        self._date = []
        self._ux = []
        self._uy = []
        self._uz = []
        self._sigux = []
        self._siguy = []
        self._siguz = []

        with open(os.path.join(prefix, filename), 'r') as file:
            reader = csv.reader(file, delimiter = ',')
            for row in reader:
                if row[0][0] == "#":
                    line = row[0]
                    match = re.search('Source File: ([A-Z0-9]{4})', line)
                    if match:
                        self._name = match.group(1)
                elif row[0] != "Datetime":
                    self._date += [isoparse(row[0])]
                    self._ux += [float(row[2])]
                    self._uy += [float(row[1])]
                    self._uz += [float(row[3])]
                    self._sigux += [float(row[5])]
                    self._siguy += [float(row[4])]
                    self._siguz += [float(row[6])]
        if self._name is not None:
            self._date = np.array(self._date, dtype = np.datetime64)
            self._ux = np.array(self._ux)
            self._uy = np.array(self._uy)
            self._uz = np.array(self._uz)
            self._sigux = np.array(self._sigux)
            self._siguy = np.array(self._siguy)
            self._siguz = np.array(self._siguz)
        else:
            raise ValueError('File was not parsed correctly.')

        from scipy.signal import detrend
        from scipy.optimize import fsolve
        from scipy.special import erf

        def find_zscore(f):
            def inner_function(zscore):
                return 1 - f - erf(zscore / np.sqrt(2.0))

            return fsolve(inner_function, 3)

        dates = self._date
        value_ux = self._ux
        value_uy = self._uy
        value_uz = self._uz
        cleared = False

        while not cleared:
            value_ux -= np.mean(value_ux)
            value_ux = detrend(value_ux)
            zscore_ux = np.abs(value_ux) / np.std(value_ux)

            value_uy -= np.mean(value_uy)
            value_uy = detrend(value_uy)
            zscore_uy = np.abs(value_uy) / np.std(value_uy)

            value_uz -= np.mean(value_uz)
            value_uz = detrend(value_uz)
            zscore_uz = np.abs(value_uz) / np.std(value_uz)

            number_of_points = value_ux.size

            f = 1/number_of_points

            max_zscore = find_zscore(f)

            i = np.where(np.logical_and(np.logical_and(zscore_ux <= max_zscore, zscore_uy < max_zscore), zscore_uz < max_zscore))

            if len(i[0]) == number_of_points:
                cleared = True
            else:
                value_ux = value_ux[i]
                value_uy = value_uy[i]
                value_uz = value_uz[i]
                dates = dates[i]

        _, indexes, _ = np.intersect1d(self._date, dates, assume_unique = False, return_indices = True)
        self._ux = self._ux[indexes]
        self._uy = self._uy[indexes]
        self._uz = self._uz[indexes]
        self._sigux = self._sigux[indexes]
        self._siguy = self._siguy[indexes]
        self._siguz = self._siguz[indexes]
        self._date = self._date[indexes]

        dates = (self._date - np.min(self._date)).astype('timedelta64[s]') / (60.0*60.0*24.0*365.0)
        G = np.ones((dates.size, 2))
        G[:,1] = dates
        try:
            invKernel = np.matmul(inv(np.matmul(G.T,G)),G.T)
            mx = np.matmul(invKernel, self._ux)
            my = np.matmul(invKernel, self._uy)
            mz = np.matmul(invKernel, self._uz)
            self._vx = mx[0]
            self._vy = my[0]
            self._vz = my[0]

            self._ux -= np.matmul(G, mx)
            self._uy -= np.matmul(G, my)
            self._uz -= np.matmul(G, mz)
        except:
            print('Problem with detrending station: {station}'.format(station = self._name))
            self._vx = 0.0
            self._vy = 0.0
            self._vz = 0.0

    def apply_filter_to_timeseries(self, filter, field_name_prefix, operation_field_prefix = None):
        if operation_field_prefix is None:
            setattr(self, field_name_prefix + '_ux', filter(self.time, self.ux))
            setattr(self, field_name_prefix + '_uy', filter(self.time, self.uy))
            setattr(self, field_name_prefix + '_uz', filter(self.time, self.uz))
            setattr(self, field_name_prefix + '_sigux', filter(self.time, self.sigux))
            setattr(self, field_name_prefix + '_siguy', filter(self.time, self.siguy))
            setattr(self, field_name_prefix + '_siguz', filter(self.time, self.siguz))
        else:
            setattr(self, field_name_prefix + '_ux', filter(self.time, getattr(self, operation_field_prefix + '_ux')))
            setattr(self, field_name_prefix + '_uy', filter(self.time, getattr(self, operation_field_prefix + '_uy')))
            setattr(self, field_name_prefix + '_uz', filter(self.time, getattr(self, operation_field_prefix + '_uz')))
            setattr(self, field_name_prefix + '_sigux', filter(self.time, getattr(self, operation_field_prefix + '_sigux')))
            setattr(self, field_name_prefix + '_siguy', filter(self.time, getattr(self, operation_field_prefix + '_siguy')))
            setattr(self, field_name_prefix + '_siguz', filter(self.time, getattr(self, operation_field_prefix + '_siguz')))

    def data_from_date(self, date, field_prefix = None):
        i = np.where(self._date == np.datetime64(date))
        if i[0].size == 0:
            return date, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        else:
            ux = float(self._ux[i] if field_prefix is None else getattr(self, field_prefix + '_ux')[i])
            uy = float(self._uy[i] if field_prefix is None else getattr(self, field_prefix + '_uy')[i])
            uz = float(self._uz[i] if field_prefix is None else getattr(self, field_prefix + '_uz')[i])
            sigux = float(self._sigux[i] if field_prefix is None else getattr(self, field_prefix + '_sigux')[i])
            siguy = float(self._siguy[i] if field_prefix is None else getattr(self, field_prefix + '_siguy')[i])
            siguz = float(self._siguz[i] if field_prefix is None else getattr(self, field_prefix + '_siguz')[i])
            return date, ux, uy, uz, sigux, siguy, siguz

    def field(self, fieldname):
        return getattr(self, fieldname)


    @property
    def name(self):
        return self._name

    @property
    def time(self):
        return self._date

    @property
    def ux(self):
        return self._ux

    @property
    def uy(self):
        return self._uy

    @property
    def uz(self):
        return self._uz

    @property
    def vx(self):
        return self._vx

    @property
    def vy(self):
        return self._vy

    @property
    def vz(self):
        return self._vz

    @property
    def sigux(self):
        return self._sigux

    @property
    def siguy(self):
        return self._siguy

    @property
    def siguz(self):
        return self._siguz


def load_all_time_series_in_path(path):
    filenames = return_all_filenames_in_path(path = path)
    return [GPSTimeSeries(filename, prefix = path) for filename in filenames]

def apply_filter_to_all_time_series(filter, time_series, field_name_prefix = 'filtered', operation_field_prefix = None):
    for ts in time_series:
        ts.apply_filter_to_timeseries(filter, field_name_prefix=field_name_prefix, operation_field_prefix=operation_field_prefix)

def gps_data_for_date_from_timeseries(date, timeseries, field_prefix = None):
    gpsdata = []
    for ts in timeseries:
        _, ux, uy, uz, sigux, siguy, siguz = ts.data_from_date(date, field_prefix=field_prefix)
        gpsdata += [GPSData(date = date, e = ux, n = uy, u = uz, se = sigux, sn = siguy, su = siguz, name = ts.name)]
    return gpsdata

def strain_from_network_for_date(strain_calculator, timeseries, date, field_prefix = None):
    gpsdata = gps_data_for_date_from_timeseries(date, timeseries, field_prefix=field_prefix)
    return strain_calculator.compute_gridded(gpsdata)

def vertical_displacements_for_date_from_timeseries(date, timeseries, network, field_prefix = None):
    lats, lons, uzs = [], [], []
    for ts in timeseries:
        _, _, _, uz, _, _, _ = ts.data_from_date(date, field_prefix=field_prefix)
        lon, lat = network.get_lat_lon_for_name(ts.name)
        lats += [lat]
        lons += [lon]
        uzs += [uz]
    return lons, lats, uzs

def check_for_nodata_at_date(date, timeseries, field_prefix = None):
    for ts in timeseries:
        name, ux, uy, uz, _, _, _ = ts.data_from_date(date, field_prefix=field_prefix)
        if np.isnan(ux) or np.isnan(uy) or np.isnan(uz):
            return False
    return True

def subsample_time_series_in_interval(time_series, min_date, max_date):
    from utils import daterange
    subset = []
    for ts in time_series:
        date_list = daterange(min_date, max_date)
        for date in date_list:
            date, ux, uy, uz, _, _, _ = ts.data_from_date(date)
            if ~np.isnan(ux) and ~np.isnan(uy) and ~np.isnan(uz):
                subset += [ts]
                break
    return subset

