import numpy as np
import csv

class EddyData():
    
    def __init__(self, filename, time_offset = np.timedelta64(-8, 'h')):
        
        time_offset = np.timedelta64(-8, 'h')
        reader = csv.reader(open(filename,'r',encoding = 'utf-8-sig'), delimiter=',')
        columns = None
        self.data_structure = {}

        for line in reader:
            if columns is None:
                columns = line[0:]
                self.data_structure = {key:[] for key in columns}
            else:
                self.data_structure[columns[0]] += [np.datetime64(line[0])+time_offset]
                for (key, data) in zip(columns[1:], line[1:]):
                    self.data_structure[key] += [np.float(data)]

        self.data_structure[columns[0]] = np.array(self.data_structure[columns[0]], dtype = np.datetime64)
        for key in columns[1:]:
            self.data_structure[key] = np.array(self.data_structure[key])

    def __apply_filtered_indexes(self, filtered_indexes):
        for key in self.data_structure:
            self.data_structure[key] = self.data_structure[key][filtered_indexes]

    def filter_by_ustar(self, ustar):
        i = np.where(np.all((self.data_structure['qc_co2_flux'] < 1, self.data_structure['ustar'] >= ustar), axis=0))
        self.__apply_filtered_indexes(i)
        
    def denoise(self, field):
        from scipy.signal import detrend
        from scipy.optimize import fsolve
        from scipy.special import erf
    
        def find_zscore(f):
            def inner_function(zscore):
                return 1 - f - erf(zscore / np.sqrt(2.0))
        
            return fsolve(inner_function, 3)
    
        value_m = self.data_structure[field]
        dates_m = self.data_structure['Datetime']
    
        cleared = False
    
        while not cleared:
            value_m = value_m - np.mean(value_m)
    
            value_m = detrend(value_m)

            zscore = np.abs(value_m) / np.std(value_m)
    
            number_of_points = value_m.size
    
            f = 1/number_of_points
        
            max_zscore = find_zscore(f)
        
            i = np.where(zscore <= max_zscore)
        
            if len(i[0]) == number_of_points:
                cleared = True
            else:
                value_m = value_m[i]
                dates_m = dates_m[i]
        
        _, indexes, _ = np.intersect1d(self.data_structure['Datetime'], dates_m, assume_unique = False, return_indices = True)
        
        self.__apply_filtered_indexes(indexes)
        
    def detrend(self, field):
        
        from scipy.signal import detrend
        return_data = self.data_structure[field]
        return_data -= np.mean(return_data)
        return_data = detrend(return_data)
        
        return return_data
        
    @property
    def time(self):
        return self.data_structure['Datetime']
    
    @property
    def co2_flux(self):
        return self.data_structure['co2_flux']
    
    