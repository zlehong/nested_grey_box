# Nested surrogate grey box energy model to run with archetypal inputs (umi tempaltes)
# Zoe Le Hong
# 11/2022
########################

# TODO:
# - sunlight calculations
# - nat vent schedule

import csv
import os

import numpy as np
import pandas as pd
import scipy.interpolate as interp

# import archetypal
# import eppy

################################
## CONSTANTS AND COEFFICIENTS ##
################################

# Coefficents of SHGC based Kmeans clustering study using Denver's climate
SHGC_DIFF_COEFF = [[17.0, 0.93], [41.0, 2.03], [55.0, 1.87], [67.0, 1.6], [78.0, 1.36], [90, 1.1]] #[[18.0, 0.83], [43.0, 1.65], [57.0, 1.58], [68.0, 1.4], [79.0, 1.23], [90, 0.98]]
SHGC_DIR_COEFF = [[34.0, 0.94], [42.0, 0.93], [51.0, 0.89], [61.0, 0.82], [70.0, 0.67], [90, 0.49]]
SHGC_79 = 0.72
SHGC_79_DIF = 1.1

# EPW indices from raw data, shifted by removing accuracy flag at index 5
IDX_LAT = 6
IDX_LNG = 7
IDX_TZ = 8
IDX_ALT = 9

IDX_CITY = 1
IDX_RE = 2
IDX_CNT = 3
IDX_EPW = 4

IDX_DBT = 5
IDX_DPT = 6
IDX_RH = 7
IDX_ATM = 8
IDX_GHR = 12
IDX_DNR = 13
IDX_DR = 14
IDX_WD = 19
IDX_WV = 20

## HEAT GAINS: POSITIVE, HEAT LOSSES: NEGATIVE

################################

################################
## LOAD IDF AND WEATHER FILE ##
################################

class Climate:
    '''
    Class for loading and parsing epw files.
    '''
    def __init__(self, epw_path):
        self.tstep = 1
        self.epw_path = epw_path
        self.epw_list = self.read_epw()
        self.epw_meta = self.parse_epw_metad()
        self.parse_epw_data()
        # self.calc_degree_days()
        # self.get_cz()
    
    def read_epw(self):
        '''
        Open epw file and store as a 
        '''
        f = open(self.epw_path)
        epw_raw = csv.reader(f)
        return list(epw_raw)
    
    def parse_epw_metad(self):
        # lon is -1*longitude of weather data
        # time zone is -15*timezone of weather data (timezone in degrees)
        self.latitude = round(float(self.epw_list[0][IDX_LAT]), 0)
        # self.lat = self.latitude
        self.longitude = round(float(self.epw_list[0][IDX_LNG])*-1.0, 0)
        # self.long = round(self.longitude*-1.0, 0)
        self.timezone = float(self.epw_list[0][IDX_TZ])
        self.tz = round(self.timezone*-15, 0)
        self.altitude = float(self.epw_list[0][IDX_ALT])
        self.city = str(self.epw_list[0][IDX_CITY])
        self.region = str(self.epw_list[0][IDX_RE])
        self.country = self.epw_list[0][IDX_CNT]
        self.epw_type = self.epw_list[0][IDX_EPW]
        self.daylight_savings()
        return self.get_clim_info()

    def parse_epw_data(self):
        ''' EPW DATA STRUCTURE
        Y, M, D, h, m, (uncertainty,) dbt, dpt, rh, pa, extrat_hzrad, extrat_dnrad,
        hz_inf_rad, ghzrad, dnrad, difhzrad, ghzlux, dnlux, difhzlux, zenith, 
        wind_dir, wind_s, sky, opaque_sky, vis, ceiling_height, wth_obs, wth_codes,
        vap_mm, opt_depth, snow, dt_snow, rain_mm, rain_t
        '''
        raw_data = np.array(self.epw_list[8:])
        self.raw_data = np.delete(raw_data, 5, 1)
        self.parse_raw(self.raw_data)
        
    
    def parse_raw(self, raw_data):
        # weekmask = [0, 1, 1, 1, 1, 1, 0]
        self.strdatetime = np.apply_along_axis('-'.join, 1, raw_data[:, 1:5])
        self.datetime = pd.date_range("2018-01-01", freq="1H", periods=8760)
        # self.datetime = np.arange('2000-01-01', '2000-12-31', dtype='datetime64[h]').astype('datetime64[m]')
        # Initialize data
        self.dbt = raw_data[:, IDX_DBT].astype(float)
        self.dpt = raw_data[:, IDX_DPT].astype(float)
        self.rh = raw_data[:, IDX_RH].astype(np.int_)
        self.atm = raw_data[:, IDX_ATM].astype(np.int_)

        self.GlobHzRad = raw_data[:, IDX_GHR].astype(np.int_)
        self.DirNormRad = raw_data[:, IDX_DNR].astype(np.int_)
        self.DiffRad = raw_data[:, IDX_DR].astype(np.int_)

        self.windDir = raw_data[:, IDX_WD].astype(float)
        self.windVel = raw_data[:, IDX_WV].astype(float)

    def daylight_savings(self):
        b = self.epw_list[4][2]
        if b == '0':
            self.DLS = False
        else:
            self.DLS = True
    
    def get_clim_info(self):
        return {"country":self.country, "region":self.region, "city":self.city, 
                "latitude":self.latitude, "longitude":self.longitude, "timezone":self.timezone,
                "daylight savings": self.DLS, "epw": self.epw_type}

    def calc_degree_days(self, HD_thresh = 18.3, CD_thresh = 10.0):
        '''
        Get heating and cooling degree days for given climate file.
        Calculating for CDD10 and HDD18 as per ASHRAE 90.1
        '''
        if HD_thresh is not float or CD_thresh is not float:
            try:
                HD_thresh = float(HD_thresh)
                CD_thresh = float(CD_thresh)
            except:
                raise TypeError("Heating and cooling degree day setpoints must be float dtype")
        
        df_days = np.reshape(self.dbt, (-1, 24))
        x = (np.max(df_days, axis=1)+np.min(df_days, axis=1))*0.5
        hdd = np.where(x < HD_thresh, HD_thresh-x, 0)
        self.HDD = round(hdd.sum(),1)
        cdd = np.where(x > CD_thresh, x-CD_thresh, 0)
        self.CDD = round(cdd.sum(), 1)

        # Turn on heating if a heating degree day and if hour is below threshold?
        # self.heatinghours = 0
        # self.coolinghours = 0

        return self.HDD, self.CDD

    def get_cz(self):
        '''
        Get ASHRAE Climate Zone according to Standard 169, Appendix A
        '''
        self.CZ=''

        if 6000<self.CDD: self.CZ='CZ0'
        elif 5000<self.CDD<=6000: self.CZ='CZ1'
        elif 3500<self.CDD<=5000: self.CZ='CZ2'
        elif self.CDD<3500 and self.HDD<=2000: self.CZ='CZ3'
        elif 2000<self.HDD<=3000 and self.CDD<3500: self.CZ='CZ4'
        elif 3000<self.HDD<=4000 and self.CDD<=3500: self.CZ='CZ5'
        elif 4000<self.HDD<=5000: self.CZ='CZ6'
        elif 5000<self.HDD<=7000: self.CZ='CZ7'
        else: self.CZ='CZ8'

        return self.CZ
    
    def set_n_tstep(self, freq):
        '''
        Set calculation frequency. Must be within 1, 2, 3, 4, 5, 6, 19, 12, 15, 20, 30
        '''
        constr = [1, 2, 3, 4, 5, 6, 19, 12, 15, 20, 30]
        if freq in constr:
            self.interpolate(freq = freq)
            self.tstep = freq
        else:
            raise ValueError('Time step must be one of: 1, 2, 3, 4, 5, 6, 19, 12, 15, 20, 30')
        return self

    def clean_data(self):
        '''
        replace missing values (99999) with interpolated value
        TODO
        '''
        pass
    
    def interpolate_wth(self, freq):
        '''
        Linear interpolation or find mean of weather data to allow for different analysis frequencies.
        '''
        num_points = int(freq*len(self.datetime))
        new_length = np.arange(num_points)
        old_length = np.arange(0, num_points, step=freq)

        # If no change, exit
        if freq == self.tstep:
            print("No change to time steps.")
            return self

        # If new timestep is larger or smaller than previous, interpolate
        else:
            # For each column, between each number add freq-1 items and interpolate
            print(f"Interpolating from {60/self.tstep} to {60/freq} min time steps.")

            interp_f = interp.interp1d(old_length, self.datetime, fill_value="extrapolate")
            interp_time = interp_f(new_length)
            self.datetime = interp_time.astype('datetime64[m]')

            interp_f2 = interp.interp1d(old_length, self.dbt, fill_value="extrapolate")
            self.dbt = interp_f2(new_length)

            interp_f2 = interp.interp1d(old_length, self.rh, fill_value="extrapolate")
            self.rh = interp_f2(new_length)

            interp_f2 = interp.interp1d(old_length, self.atm, fill_value="extrapolate")
            self.atm = interp_f2(new_length)

            interp_f2 = interp.interp1d(old_length, self.GlobHzRad, fill_value="extrapolate")
            self.GlobHzRad = interp_f2(new_length)

            interp_f2 = interp.interp1d(old_length, self.DirNormRad, fill_value="extrapolate")
            self.DirNormRad = interp_f2(new_length)

            interp_f2 = interp.interp1d(old_length, self.DiffRad, fill_value="extrapolate")
            self.DiffRad = interp_f2(new_length)

            interp_f2 = interp.interp1d(old_length, self.windDir, fill_value="extrapolate")
            self.windDir = interp_f2(new_length)

            interp_f2 = interp.interp1d(old_length, self.windVel, fill_value="extrapolate")
            self.windVel = interp_f2(new_length)
            
            return self
    
    @classmethod
    def nat_vent_scheduler(cls):
        '''
        Determine natural ventilation "opportunity schedule." 
        Flags for natural ventilation on/off. Returns infiltration heat gains/losses.
        '''
        pass

    @classmethod
    def sun_path(cls):
        '''
        Determine solar angle
        Flags for natural ventilation on/off. Returns infiltration heat gains/losses.
        '''
        pass

class StochasticScheduler():
    '''
    To work with umi/archetype schedules and IDFs
    Change stochastically, interpolate timestep to match chosen. (Probably not needed for ML training.)
    Note, archetypal does have a random scheduler.
    '''
    pass


if __name__ == "__main__":
    path = os.getcwd()
    # epw_path = os.path.join('GreyBoxTestData', 'epw', 'USA_MN_Minneapolis-St.Paul.Intl.AP.726580_TMYx.epw')
    epw_path = os.path.join('GreyBoxTestData', 'epw', 'USA_AZ_Phoenix-Sky.Harbor.Intl.AP.722780_TMYx.epw')
    phx_wth = Climate(epw_path)
    print(phx_wth.epw_meta)
    # hdd, cdd = phx_wth.calc_degree_days()
    # print(phx_wth.CZ)
    # print(len(phx_wth.datetime))
    print(phx_wth.dbt)
    phx_wth.set_n_tstep(4)
    phx_wth.set_n_tstep(2)
    # print(phx_wth.tstep)
    print(phx_wth.dbt)

## INSTANTANEOUS GAIN CALCULATIONS ##

# class InstGains():
#     '''
#     Class for calculating instantaneous gains. 
#     '''
#     def __init__(self, epw_path):
#         pass


# def to_dict(self):
#     '''this is what I want to do'''
#     d={}
#     for v in my_python_obj.attributes:
#         d[v] = self.v
#     return d