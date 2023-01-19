import numpy as np

# from nested-grey-box import Climate
from climate import Climate

# from archetypal import UmiTemplateLibrary
# from pyumi.shoeboxer import ShoeBox

################################
## SHOEBOX GEOMETRY SUMMARY ##
################################
# Use archetypal and umi shoebox geometry definitions to organize test shoebox data
# Calculate geometric inputs (overall R-vals)
# Create idfs for each test box
"""
(total) facade area
core area
roof area (only exposed roof area?)
base area (total floor of adiabatic + ground)
adiabatic area
floor 2 facade - input from shoeboxer
core 2 perim
roof 2 floor
ground 2 floor
"""


class NGBShoebox:
    """
    Class for test shoebox data and creation of .idfs for E+ testing.
    Creates templates for test data distributions/ranges.
    Intitializes idfs from templates with archetypal.
    """

    def __init__(
        self,
        epw_path,
        # Shoebox geometry info
        width=3.0,
        height=3.0,
        length=10.5,
        perim_2_a=0.5,
        rf_2_a=0.1,
        gnd_2_a=0.1,
        weight=1,
        orientation=180,
        wwr=0.4,
    ):
        self.epw_path = epw_path
        self.climate = Climate(self.epw_path)

        # Geometric data
        self.width = width  # meters
        self.height = height
        self.length = length
        self.area = self.width * self.height
        self.perim_ratio = perim_2_a
        self.roof_ratio = rf_2_a
        self.gnd_ratio = gnd_2_a
        self.orientation = orientation  # Orientation, N = 0, E = 90, S = 180, W = 270
        self.wwr = wwr  # wwr for exposed orientation
        # self.shading

        # Shoebox weight information
        self.weight = weight

    def energy_setup(self):
        """
        Set up energy model data
        (in the sb class as some buildings may have muliple uses in future)
        """
        self._get_energy_dat()
        self._get_envelope_dat()
        self._get_schedule_dat()

    def _get_energy_dat(self):
        """
        Set up baseline energy information based on geometry, climate, archetype
        """
        self.LPD = 4  # W/m2
        self.EPD = 6  # W/m2
        self.PPD = 0.01  # people/sq.m
        self.COP_heat = 1
        self.COP_cool = 3
        self.setpoint_l = 20
        self.setpoint_u = 24
        self.setpoint_nv = 18
        self.infil = 18  # ACH50
        self.NV_ach_thresh = 10  # NVach has to be >1 for NV to be used
        # Local electric/gas CO2 and cost data
        self.CO2_elec = 0.275
        self.cost_elec = 0.26
        self.CO2_gas = 0.18
        self.cost_gas = 0.04
        # Locational data - IS THIS NECESSARY?
        self.lat = self.climate.latitude
        self.long = self.climate.longitude
        self.timezone = self.climate.timezone
        self.grnd_temp = 10  # constant or list for each month

    def _get_envelope_dat(self):
        """
        Get envelope data from idf
        (U_wall, U_roof, U_ground, U_glazing, SHGC, thermal mass factors, exterior convection convection)
        """
        self.SHGC = 0.65
        # Coefficents of SHGC based Kmeans clustering study using Denver's climate
        self.SHGCcoeff_dif = [
            [17.0, 0.93],
            [41.0, 2.03],
            [55.0, 1.87],
            [67.0, 1.6],
            [78.0, 1.36],
            [90, 1.1],
        ]  # [[18.0, 0.83], [43.0, 1.65], [57.0, 1.58], [68.0, 1.4], [79.0, 1.23], [90, 0.98]]
        self.SHGCcoeff_dir = [
            [34.0, 0.94],
            [42.0, 0.93],
            [51.0, 0.89],
            [61.0, 0.82],
            [70.0, 0.67],
            [90, 0.49],
        ]

    def _get_schedule_dat(self):
        """
        Get schedule time-series data from idf (TODO)
        """
        pass

    def set_orientation(self, angle):
        angle = angle / 360
        if angle > 315 and angle <= 45:
            self.orientation = 0
        elif angle > 45 and angle <= 135:
            self.orientation = 90
        elif angle < 135 and angle <= 225:
            self.orientation = 180
        elif angle > 225 and angle <= 315:
            self.orientation = 270

    # def templates_iter(self, template):
    #     yield next(iter(template)).BuildingTemplates

    def _gen_test_schedules(self):
        """
        Hourly schedules from Alpha Arsano's climabox for comparison
        For an office building.
        """
        ### Occupancy Schedule
        sch_occ = []
        for i in range(1, 7):
            sch_occ.append((i, 0))
        sch_occ.append((7, 0.1))
        sch_occ.append((8, 0.2))
        for i in range(9, 13):
            sch_occ.append((i, 0.95))
        sch_occ.append((13, 0.5))
        for i in range(14, 18):
            sch_occ.append((i, 0.95))
        sch_occ.append((18, 0.3))
        sch_occ.append((19, 0.1))
        sch_occ.append((20, 0.1))
        for i in range(21, 25):
            sch_occ.append((i, 0.05))

        sch_occOff_Sat = []
        for i in range(1, 8):
            sch_occOff_Sat.append((i, 0))
        sch_occOff_Sat.append((8, 0.1))
        for i in range(9, 13):
            sch_occOff_Sat.append((i, 0.3))
        for i in range(13, 18):
            sch_occOff_Sat.append((i, 0.1))
        for i in range(18, 25):
            sch_occOff_Sat.append((i, 0.0))

        # No occupancy, off schedule
        sch_Off_Sun = []
        for i in range(1, 25):
            sch_Off_Sun.append((i, 0))

        self.occSch = sch_Off_Sun + sch_occ * 5 + sch_occOff_Sat
        self.occSch = np.array(self.occSch * 52 + sch_Off_Sun, dtype=np.float16)

        ### Light schedule
        sch_lgtOff = []
        for i in range(1, 7):
            sch_lgtOff.append((i, 0.05))
        sch_lgtOff.append((7, 0.1))
        sch_lgtOff.append((8, 0.3))
        for i in range(9, 17):
            sch_lgtOff.append((i, 0.90))
        sch_lgtOff.append((17, 0.5))
        sch_lgtOff.append((18, 0.5))
        for i in range(19, 21):
            sch_lgtOff.append((i, 0.30))
        for i in range(21, 23):
            sch_lgtOff.append((i, 0.20))
        sch_lgtOff.append((23, 0.1))
        sch_lgtOff.append((24, 0.05))

        sch_lgtOff_Sat = []
        for i in range(1, 7):
            sch_lgtOff_Sat.append((i, 0.05))
        for i in range(7, 9):
            sch_lgtOff_Sat.append((i, 0.1))
        for i in range(9, 13):
            sch_lgtOff_Sat.append((i, 0.3))
        for i in range(13, 18):
            sch_lgtOff_Sat.append((i, 0.15))
        for i in range(18, 25):
            sch_lgtOff_Sat.append((i, 0.05))

        sch_lgtOff_Sun = []
        for i in range(1, 25):
            sch_lgtOff_Sun.append((i, 0.05))

        self.lgtSch = sch_lgtOff_Sun + sch_lgtOff * 5 + sch_lgtOff_Sat
        self.lgtSch = np.array(self.lgtSch * 52 + sch_lgtOff_Sun, dtype=np.float16)

        # Equipment schedule
        sch_eqpOff = []
        for i in range(1, 8):
            sch_eqpOff.append((i, 0))
        sch_eqpOff.append((8, 0.4))
        for i in range(9, 13):
            sch_eqpOff.append((i, 0.90))
        sch_eqpOff.append((13, 0.8))
        for i in range(14, 18):
            sch_eqpOff.append((i, 0.90))
        sch_eqpOff.append((18, 0.5))
        for i in range(19, 25):
            sch_eqpOff.append((i, 0.4))
        sch_l = [x[1] for x in sch_eqpOff]
        # print( sch_l)

        sch_eqpOff_Sat = []
        for i in range(1, 7):
            sch_eqpOff_Sat.append((i, 0.3))
        for i in range(7, 8):
            sch_eqpOff_Sat.append((i, 0.4))
        for i in range(8, 13):
            sch_eqpOff_Sat.append((i, 0.5))
        for i in range(13, 18):
            sch_eqpOff_Sat.append((i, 0.35))
        for i in range(18, 25):
            sch_eqpOff_Sat.append((i, 0.30))

        sch_eqpOff_Sun = []
        for i in range(1, 25):
            sch_eqpOff_Sun.append((i, 0.3))

        self.eqpSch = sch_eqpOff_Sun + sch_eqpOff * 5 + sch_eqpOff_Sat
        self.eqpSch = np.array(self.eqpSch * 52 + sch_eqpOff_Sun, dtype=np.float16)

        # sch_L_E = []
        # for i in range(1,8): sch_L_E.append((i,0.1))
        # for i in range(8,12): sch_L_E.append((i,0.2))
        # for i in range(12,19): sch_L_E.append((i,0.5))
        # for i in range(19,25): sch_L_E.append((i,0.1))

        # # NAT VENT SCHEDULES
        # sch_nightV = []
        # for i in range(1,6): sch_nightV.append((i,1))
        # for i in range(6,22): sch_nightV.append((i,0))
        # for i in range(22,25): sch_nightV.append((i,1))

        # sch_preCool = []
        # for i in range(1,6): sch_preCool.append((i,0))
        # for i in range(6,9): sch_preCool.append((i,1))
        # for i in range(9,18): sch_preCool.append((i,0))
        # for i in range(18,22): sch_preCool.append((i,1))
        # for i in range(22,25): sch_preCool.append((i,0))

    @classmethod
    def from_shoeboxer(cls, eplus_path, h=3, w=3):
        """
        Make class instance from shoeboxer idf outputs
        """
        pass

    @classmethod
    def save_idf_from_template(
        cls,
        name,
        template,
        zone_data,
        wwr_map,
        path="C:/Users/zoele/pyumi/pyumi/TrainingData/idfs/",
    ):
        """
        zone_data = [{
                "name":np.nan,
                "coordinates": [(10, 0), (10, 10), (0, 10), (0, 0)],
                "height": 3,
                "num_stories": 1,
                "zoning": "core/perim",
                "perim_depth": 3}]
        wwr_map= { # N is 0, E is 90
                0: 0,
                90: 0,
                180: 30,
                270: 0}
        """
        pass

        # sb = ShoeBox.from_template(
        #     template,
        #     zones_data=zone_data,
        #     wwr_map=wwr_map,
        #     name=name
        # )

        # Ensure calculation of radiation on floor
        # Ensure 15 min time step

        # Eppy: saveas(filename, lineendings='default', encoding='latin-1')
        # sb.saveas(path+name)

        # return sb
