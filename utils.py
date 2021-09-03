from os import listdir
from os.path import isfile, join
from datetime import timedelta

def return_all_filenames_in_path(path = '.'):
    return [f for f in listdir(path) if isfile(join(path, f))]

def daterange(start_date, end_date, days = 1):
    for n in range(int((end_date - start_date).total_seconds() / (60.0 * 60.0 * 24.0 * days))):
        yield start_date + timedelta(seconds=days*n*60.0*60.0*24.0)