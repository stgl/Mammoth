import numpy as np

def gaussian_filter(wavelength_in_days):

    def inner_function(times, data):

        time_decimal_day = (times - np.min(times)).astype('timedelta64[s]').astype(np.float64) / (60.0*60.0*24.0)

        def calculate_values(t0):
            tn = (time_decimal_day - t0) / wavelength_in_days
            i = np.where(np.logical_and(tn > -5, tn < 5))
            weights = np.exp(-0.5*np.power(tn[i],2))
            weights /= np.sum(weights)
            return np.sum(data[i]*weights)

        return np.array([calculate_values(t) for t in time_decimal_day])

    return inner_function

def reference_to_date(date):

    def inner_function(times, data):
        i = np.where(np.datetime64(date) == times)
        return data - data[i]

    return inner_function
