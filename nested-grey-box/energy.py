import numpy as np
import pandas as pd

import math

# from NestedGreyBox.climate import Climate
# from NestedGreyBox.geometry import NGBShoebox

from pvlib import irradiance
# from pvlib import clearsky, atmosphere, solarposition
from pvlib.location import Location
from pvlib.solarposition import declination_spencer71

################################
## OUTPUTS TO MATCH UMI IDFS ##
################################
'''
heating
cooling
----
lighting
people_gain
solar_gain
infil_gain
infil_lossclimabox_ng
vent_loss
vent_gain
net_vent_gain
nat_vent_loss
electric_equip
gas_equip
hot_water
'''

class RadFacade:
    '''
    Adapted from Alpha Arsano's ClimaBox and Sam D. code to calculate facade radiation for exterior surface.
    '''

    def __init__(self, clim_obj, solar_position_method='nrel_numpy'):
        '''
        Init with weather in pd database format with ['ghi', 'dni', 'dhi']
        Location data (city, lat, long, time zone)
        solar_position_method (str): default 'nrel_numpy'.
                'nrel_numpy' uses an implementation of the NREL SPA algorithm.
                'nrel_numba' uses an implementation of the NREL SPA algorithm, but also compiles the code first.
                'pyephem' uses the PyEphem package.
                'ephemeris' uses the pvlib ephemeris code. 'nrel_c' uses the NREL SPA C code.
        '''
        # TODO: altitude to be inherited from shoeboxer (shoebox altitude, not weather file)

        self.clim_obj = clim_obj

        pv_tz = f"Etc/GMT{int(self.clim_obj.timezone):+d}"

        self.location = Location(
            self.clim_obj.latitude, self.clim_obj.longitude, name=self.clim_obj.city, altitude=self.clim_obj.altitude, tz=pv_tz
        )
        self.solar_position_method = solar_position_method

        weather = pd.DataFrame(columns = ["ghi", "dni", "dhi"])
        weather['ghi'] = self.clim_obj.GlobHzRad
        weather['dni'] = self.clim_obj.DirNormRad
        weather['dhi'] = self.clim_obj.DiffRad
        idx = pd.DatetimeIndex(self.clim_obj.datetime).tz_localize(self.location.tz)
        weather = weather.set_index(idx)
        self.weather = weather

        self.declination = np.vectorize(declination_spencer71)(
            self.weather.index.dayofyear
        )
        self.solar_position = self.location.get_solarposition(
            self.weather.index, method=self.solar_position_method
        )
        self.clearsky = self.location.get_clearsky(self.weather.index)

        # clim_obj.datetime - this is index for input to radfacade
        # repeat timezone in db input for rad facade

    @property
    def azimuth(self):
        return self.solar_position["azimuth"].values

    @property
    def location(self):
        """Get or set the location."""
        return self._location

    @location.setter
    def location(self, value):
        assert isinstance(value, Location)
        self._location = value

    @property
    def weather(self) -> pd.DataFrame:
        """Get or set the weather DataFrame.
        Note: columns are ``ghi, dni, dhi``.
        """
        return self._weather

    @weather.setter
    def weather(self, value):
        assert isinstance(value, pd.DataFrame), "weather must be a DataFrame."
        assert sorted(value.columns.tolist()) == sorted(["ghi", "dni", "dhi"]), (
            f"columns of weather must be '['ghi', 'dni', 'dhi']', "
            f"not '{value.columns.tolist()}'"
        )
        assert isinstance(
            value.index, pd.DatetimeIndex
        ), f"weather.index must be a '{pd.DatetimeIndex}', not a '{type(value.index)}'"
        assert (
            value.index.tzinfo is not None
        ), "weather.index must be localized. use 'pandas.DatetimeIndex.tz_localize'."
        self._weather = value

    @property
    def solar_position_method(self):
        """Get or set the solar position method."""
        return self._solar_position_method

    @solar_position_method.setter
    def solar_position_method(self, value):
        assert value in [
            "nrel_numpy",
            "nrel_numba",
            "pyephem",
            "ephemeris",
            "nrel_c",
        ], (
            "Input value error for '{value}'. solar_position_method must be one of ("
            "'nrel_numpy', 'nrel_numba', 'pyephem', 'ephemeris', 'nrel_c')"
        )
        self._solar_position_method = value

    @property
    def solar_position(self):
        """Calculate the solar zenith, azimuth, etc. at this location."""
        return self._solar_position

    @solar_position.setter
    def solar_position(self, value):
        self._solar_position = value

    def get_irradiance_for_surface(self, tilt, surface_azimuth, times=None):
        """Calculate clear-sky GHI and transpose to plane of array Define a function
        so that we can re-use the sequence of operations with different locations.
        Args:
            tilt (float):
            surface_azimuth (float): surface azimuth from north.
        Returns:
            pd.DataFrame: total_irrad - Contains keys/columns 'poa_global',
                'poa_direct', 'poa_diffuse', 'poa_sky_diffuse', 'poa_ground_diffuse'.
        """
        if times is None:
            times = self.weather.index

        # Get weather data using values for GHI, DNI, and DHI
        # weather = self.weather.loc[times, :]
        weather = self.weather.loc[times, :]
        # Get solar azimuth and zenith to pass to the transposition function
        solar_position = self.solar_position.loc[times, :]
        # Use the get_total_irradiance function to transpose the GHI to POA
        poa_irradiance = irradiance.get_total_irradiance(
            surface_tilt=tilt,
            surface_azimuth=surface_azimuth,
            dni=weather["dni"],
            ghi=weather["ghi"],
            dhi=weather["dhi"],
            solar_zenith=solar_position["zenith"],
            solar_azimuth=solar_position["azimuth"],
        )
        # Return DataFrame
        
        return poa_irradiance # Units are W/m2

class ShoeboxEnergy:
    '''
    Class for calculating energy use of entire building (can be multiple shoeboxes)
    '''
    def __init__(self, geom_obj, SHGCbool=True):
        self.shoebox = geom_obj
        self.climate = self.shoebox.climate
        self.orientation = self.shoebox.orientation
        self.SHGCbool = SHGCbool
    
    def _get_tsd(self):
        '''
        Set up numpy arrays for energy calculations
        '''
        # Climate data
        self.dbt = self.climate.dbt
        # Schedule data
    
    def calc_solar_gains(self):
        # TODO: make this work for windows in all directions
        self._rad_obj = RadFacade(self.climate)
        self._dirrad_facade, self._difrad_facade = (self._rad_obj.get_irradiance_for_surface(tilt = 90, surface_azimuth = self.orientation).loc[:, ["poa_direct", "poa_diffuse"]].values.T)
        self.dir_sol_gain = self.shoebox.width * self.shoebox.height * self.shoebox.wwr * self._dirrad_facade # Units: W
        self.dif_sol_gain = self.shoebox.width * self.shoebox.height * self.shoebox.wwr * self._difrad_facade
        
        # Adjust for SHGC based on hourly angle of incidence (AOI)
        # Based on Alpha Arsano ClimaBox
        shgc__79_dir = 0.72
        shgc__79_dif = 0.85

        # AOI relative to ground plane - each hour
        AOI = [irradiance.aoi(0, self.orientation, z, a) for z, a in zip(self._rad_obj.solar_position["zenith"], self._rad_obj.solar_position["azimuth"])]
        AOI = np.array(AOI)
        AOI[AOI < -90] = 0
        AOI[AOI > 90] = 0
        # QUESTION: WHY does rf.dirradS[i]>0.5???
        AOI[np.where(self._dirrad_facade < 0.5)] = 0
        # print(len(AOI))

        shgc_tempDif = shgc_tempDir = np.zeros(len(AOI))

        for a in range (len(self.shoebox.SHGCcoeff_dir)):
            shgc_tempDir[np.where(AOI < self.shoebox.SHGCcoeff_dir[a][0])] = self.shoebox.SHGCcoeff_dir[a][1]*shgc__79_dir
            # if AOI < self.shoebox.SHGCcoeff_dir[a][0]: 
            #     shgc_tempDir = self.shoebox.SHGCcoeff_dir[a][1]*shgc__79_dir
            #     break

        for a in range (len(self.shoebox.SHGCcoeff_dif)):
            shgc_tempDif[np.where(AOI < self.shoebox.SHGCcoeff_dif[a][0])] = self.shoebox.SHGCcoeff_dif[a][1]*shgc__79_dif
            # if AOI < self.shoebox.SHGCcoeff_dif[a][0]: 
            #     shgc_tempDif = self.shoebox.SHGCcoeff_dif[a][1]*shgc__79_dif
            #     break

        # TODO: Consider shading based on depth

        # If varSHGCbool==0, static SHGC of selected glazing is used. 
        # Else if varSHGCbool is 1, then variable SHGC based on solar angle is used.
        temp_shading = 0.8
        if self.SHGCbool == False:
            solarGain = self.dir_sol_gain*self.shoebox.SHGC*temp_shading + self.dif_sol_gain*self.shoebox.SHGC*temp_shading
        else: 
            solarGain = self.dir_sol_gain*shgc_tempDir*temp_shading + self.dif_sol_gain*shgc_tempDif*temp_shading
        self.SolarGain = solarGain
        # self.SolarGain = solarGain/(self.shoebox.width * self.shoebox.height)
        # self.SolarGain = [round(x/1000, 2) for x in solarGain]
        return self.SolarGain # units are W/sq.m.

    def calc_internal_gains(self):
        # W/m2 x schedule % on -- units are W/sq.m.
        self.LightGain_hr = self.shoebox.LPD*self.shoebox.lgtSch[:,1] #*LPD_floorPercentage
        self.EquipGain_hr = self.shoebox.EPD*self.shoebox.eqpSch[:,1]
        self.PplGain_hr = self.shoebox.PPD*self.shoebox.occSch[:,1]*self.shoebox.width*self.shoebox.length