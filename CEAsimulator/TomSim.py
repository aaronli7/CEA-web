import csv, json, time
import numpy as np
import os, configparser, math
from .tools.utils import *

# constant definition
PI = 3.141592654
LAT = 52.0
# self.IAGE = [0] * 75


class TomSim:
    """
    This class simulates the tomato production.

    params:
    - co2: CO2 in the greenhouse.
    - temperature: temperature in the greenhouse.
    - fruit_per_truss: number of fruits in one truss.
    - lon: longitude of the site.
    - lat: latitude of the site.
    - start_date: the julian day of the simulation start date.
    - end_date: the julian day of the simulation end date.
    """

    def __init__(self, co2, temperature, fruit_per_truss, lon, lat, start_date, end_date, debug=False):
        # user define
        self.temperature = temperature
        self.co2 = co2
        self.dbg = debug
        self.fruit_per_truss = fruit_per_truss
        self.lon = lon
        self.lat = lat
        self._start_time = start_date
        self._finish_time = end_date

        self.IAGE = [0] * 75
        self._hourly_radiation_table = np.zeros((730, 24))
        self._current_day_radiation = np.zeros(24)
        self._hourly_temperature_table = np.full((730, 24), self.temperature)
        self._current_day_temperature = np.zeros(24)
        self._hourly_co2_table = np.full((730, 24), self.co2)
        self._current_day_co2 = np.zeros(24)
        self._daily_radiation_sum_table = np.zeros(730)
        self._daily_average_co2_table = np.zeros(730)
        self._daily_average_temperature_table = np.zeros(730)
        self._fraction_dry_matter_leaves_table = np.zeros(730)
        self._fraction_dry_matter_stem_table = np.zeros(730)
        self._fraction_dry_matter_organs_table = np.zeros(730)
        self._fraction_dry_matter_root_table = np.zeros(730)
        self._LAI_table = np.zeros(730)
        self._SLA_table = np.zeros(730)
        self._leaves_removal_table = np.zeros(730)
        self._stems_removal_table = np.zeros(730)
        self._organs_removal_table = np.zeros(730)
        self.growth_rate_table = np.zeros(730)
        self._PDVU = np.zeros(730)  # qili: don't understand the meaning
        self.fruits_per_truss_table = np.zeros(75)
        self.veg_dry_weight_unit = np.zeros(75)
        self.truss_weight_table = np.zeros(75)
        self.if_truss_harvested = np.zeros(75)
        self.veg_growth = np.zeros(75)
        self.sink = np.zeros(75)
        self.veg_sink = np.zeros(75)
        self.itg_veg = np.zeros(75)
        self.dev_stages_per_truss = np.zeros(75)
        self.growth_truss = np.zeros(75)
        self.temp_array = np.zeros(200)
        self.avtemp_truss = np.zeros(75)
        self.veg_dev = np.zeros(75)

        self._reference_year = True
        self._calc_LAI_from_SLA = True
        self._calc_SLA_standard = True
        self._fixed_crop = False
        self._manual_leaf_picking = False
        self._side_shoot = False
        self._initiation = False

        self.truss_pos = 100
        self.veg_pos = 100
        self.buffer = 0
        self.wtruss_615 = 0
        self.num_veg_units = 0

    def read_radiation(self):

        # reading the radiation from API
        df = query_ghi(lon=self.lon, lat=self.lat)
        for _, row in df.iterrows():
            [num_day, num_hours, radiation] = [
                int(row['Day']),
                int(row['Hour']),
                float(row['GHI']),
            ]
            if num_hours == 0:
                rad_sum = 0
                co2_sum = 0
                temperature_sum = 0
                count = 0
            self._hourly_radiation_table[num_day - 1, num_hours] = radiation
            rad_sum += radiation
            co2_sum += self.co2
            temperature_sum += self.temperature
            count += 1  # qili: what does this use for? it's always 24
            if num_hours == 23:
                rad_sum *= 0.36  # qili: magic number?
                self._daily_radiation_sum_table[num_day - 1] = rad_sum
                self._daily_average_temperature_table[num_day - 1] = (
                    temperature_sum / 24
                )
                self._daily_average_co2_table[num_day - 1] = co2_sum / count

        if self.dbg:
            print("Debug msg:")
            print((self._hourly_radiation_table.shape))
        return

    def read_config(self, filename):
        # config_file = os.path.abspath(os.path.join(os.getcwd(), filename))
        # print(config_file)
        config_file = filename
        cf = configparser.ConfigParser()
        cf.read(config_file)
        print("Config file path:", config_file) if self.dbg else None

        self._exti_coef_diff_light = json.loads(
            cf.get("simulation", "exti_coef_diff_light")
        )
        self._scat_coef_leaves = json.loads(cf.get("simulation", "scat_coef_leaves"))
        self._maint_respiration_coef_leaves = json.loads(
            cf.get("simulation", "maint_respiration_coef_leaves")
        )
        self._maint_respiration_coef_stems = json.loads(
            cf.get("simulation", "maint_respiration_coef_stems")
        )
        self._maint_respiration_coef_organs = json.loads(
            cf.get("simulation", "maint_respiration_coef_organs")
        )
        self._maint_respiration_coef_roots = json.loads(
            cf.get("simulation", "maint_respiration_coef_roots")
        )
        self._Q10_value = json.loads(cf.get("simulation", "Q10_value"))
        self._ref_temp = json.loads(
            cf.get("simulation", "ref_temperature_for_maint_respiration_coef")
        )
        self._assimilate_req_leaves = json.loads(
            cf.get("simulation", "assimilate_req_leaves")
        )
        self._assimilate_req_stem = json.loads(
            cf.get("simulation", "assimilate_req_stem")
        )
        self._assimilate_req_organs = json.loads(
            cf.get("simulation", "assimilate_req_organs")
        )
        self._assimilate_req_roots = json.loads(
            cf.get("simulation", "assimilate_req_roots")
        )
        self._init_dry_weight_leaves = json.loads(
            cf.get("simulation", "init_dry_weight_leaves")
        )
        self._init_dry_weight_stem = json.loads(
            cf.get("simulation", "init_dry_weight_stem")
        )
        self._init_dry_weight_orgrans = json.loads(
            cf.get("simulation", "init_dry_weight_orgrans")
        )
        self._init_dry_weight_roots = json.loads(
            cf.get("simulation", "init_dry_weight_roots")
        )
        self._plant_density = json.loads(cf.get("simulation", "plant_density"))
        # self._start_time = json.loads(cf.get("simulation", "start_time"))
        # self._finish_time = json.loads(cf.get("simulation", "finish_time"))
        self._sim_time_step = json.loads(cf.get("simulation", "sim_time_step"))
        self._output_time_step = json.loads(cf.get("simulation", "output_time_step"))
        self._frac_dry_weight_leaves = json.loads(
            cf.get("simulation", "frac_dry_weight_leaves")
        )
        self.LAI_table_size = json.loads(cf.get("simulation", "LAI_table_size"))
        self.SLA_table_size = json.loads(cf.get("simulation", "SLA_table_size"))
        self.PDVU_table_size = json.loads(cf.get("simulation", "PDVU_table_size"))
        # qili: some magic number
        (
            self._fraction_dry_matter_leaves_table[0],
            self._fraction_dry_matter_leaves_table[1],
            self._fraction_dry_matter_leaves_table[2],
            self._fraction_dry_matter_leaves_table[3],
        ) = (1, 0, 365, 0)
        (
            self._fraction_dry_matter_stem_table[0],
            self._fraction_dry_matter_stem_table[1],
            self._fraction_dry_matter_stem_table[2],
            self._fraction_dry_matter_stem_table[3],
        ) = (1, 0, 365, 0)
        (
            self._fraction_dry_matter_organs_table[0],
            self._fraction_dry_matter_organs_table[1],
            self._fraction_dry_matter_organs_table[2],
            self._fraction_dry_matter_organs_table[3],
        ) = (1, 0, 365, 0)
        (
            self._fraction_dry_matter_root_table[0],
            self._fraction_dry_matter_root_table[1],
            self._fraction_dry_matter_root_table[2],
            self._fraction_dry_matter_root_table[3],
        ) = (1, 0, 365, 0)
        (
            self._LAI_table[0],
            self._LAI_table[1],
            self._LAI_table[2],
            self._LAI_table[3],
        ) = (1, 0, 365, 0)
        (
            self._SLA_table[0],
            self._SLA_table[1],
            self._SLA_table[2],
            self._SLA_table[3],
        ) = (1, 0, 365, 0)
        (
            self._leaves_removal_table[0],
            self._leaves_removal_table[1],
            self._leaves_removal_table[2],
            self._leaves_removal_table[3],
        ) = (1, 0, 365, 0)
        (
            self._stems_removal_table[0],
            self._stems_removal_table[1],
            self._stems_removal_table[2],
            self._stems_removal_table[3],
        ) = (1, 0, 365, 0)
        (
            self._organs_removal_table[0],
            self._organs_removal_table[1],
            self._organs_removal_table[2],
            self._organs_removal_table[3],
        ) = (1, 0, 365, 0)
        (
            self._PDVU[0],
            self._PDVU[1],
            self._PDVU[2],
            self._PDVU[3],
        ) = (1, 0, 365, 0)

        self._initial_truss = json.loads(cf.get("simulation", "initial_truss"))
        self._start_day_of_truss_growth = json.loads(
            cf.get("simulation", "start_day_of_truss_growth")
        )
        self._sink_factor = json.loads(cf.get("simulation", "sink_factor"))

    def read_greenhouse_config(self, filename):
        # config_file = os.path.abspath(os.path.join(os.getcwd(), filename))
        config_file = filename
        cf = configparser.ConfigParser()
        cf.read(config_file)
        self.diffuse_rad_transmission = json.loads(
            cf.get("greenhouse", "diffuse_rad_transmission")
        )
        self.deviation_orient = json.loads(cf.get("greenhouse", "deviation_orient"))
        self.num_layer = json.loads(cf.get("greenhouse", "num_layer"))
        self.elevation_layer = json.loads(cf.get("greenhouse", "elevation_layer"))
        self.num_entries = np.asarray(json.loads(cf.get("greenhouse", "num_entries")))
        self.azimuth = np.asarray(json.loads(cf.get("greenhouse", "azimuth")))
        self.construction_transmission = np.asarray(
            json.loads(cf.get("greenhouse", "construction_transmission"))
        )
        self.glass_transmission = np.asarray(
            json.loads(cf.get("greenhouse", "glass_transmission"))
        )
        print(
            "glass transmission matrix:\n",
            self.glass_transmission,
        ) if self.dbg else None

    def pre_simulation(self):
        self.time = self._start_time
        self.icount = 1
        self.icount1 = 1

        # correction on number of fruit per truss
        self.cor_fruit_per_truss = 1.0
        # correction on leaf picking
        self.cor_leaf_picking = 1.0

        for i in range(50):
            self.fruits_per_truss_table[i] = (
                self.cor_fruit_per_truss * self.fruit_per_truss
            )

        self.total_dry_weight = (
            self._init_dry_weight_leaves
            + self._init_dry_weight_stem
            + self._init_dry_weight_orgrans
            + self._init_dry_weight_roots
        )
        self.dry_weight_leaves_existing = self._init_dry_weight_leaves
        self.dry_weight_stems_existing = self._init_dry_weight_stem
        self.dry_weight_organs_existing = self._init_dry_weight_orgrans
        self.truss_weight_table[0] = self._init_dry_weight_orgrans

        self.veg_dry_weight_unit[0] = (
            self._init_dry_weight_leaves + self._init_dry_weight_stem
        ) * 0.7
        self.veg_dry_weight_unit[1] = (
            self._init_dry_weight_leaves + self._init_dry_weight_stem
        ) * 0.2
        self.veg_dry_weight_unit[2] = (
            self._init_dry_weight_leaves + self._init_dry_weight_stem
        ) * 0.1
        self.RGRSOM = 0  # qili: I don't know the meaning of it.

        # qili: fixed parameters, I don't know the meaning of them.
        self.A = 1.0
        self.B = 33.0

        # correction on greenhouse transmissivity
        self.cor_transmissivity = 1.13
        # correction on maintenance respiration
        self.cor_maintenance_respiration = 1.0
        # correction on growth respiration
        self.cor_growth_respiration = 1.0
        # correction on vegetative sink
        self.cor_veg_sink = 1.0
        # correction on irradiance
        self.cor_irradiance = 1.0
        # correction on temperature
        self.cor_temperature = 1.0
        # correction on CO2
        self.cor_co2 = 1.0
        # correction on plant density
        self.cor_plant_density = 1.0
        # correction on SLA (SLA-function normalized)
        self.cor_SLA = 1.0
        # correction on flowering rate
        self.cor_flower_rate = 1.0
        # correction on fruit development stage
        self.cor_fruit_dev = 1.0
        # ground reflection coefficient
        self.ground_ref_coef = 0.5

        if self._side_shoot:
            self._flower_shoot_date = 100
            self.new_plant_density = 3.125
            self._plant_density = new_plant_density
        self._plant_density *= self.cor_plant_density

    def daytime_radiation_character(self, day, lat_rad):
        rad_factor = 0.017453292
        if lat_rad > 67 * rad_factor:
            quit("ERROR in ASTROG: LAT > 67")
        elif lat_rad < -67 * rad_factor:
            quit("ERROR in ASTROG: LAT < - 67")

        self.sun_declination = -math.asin(
            math.sin(23.45 * rad_factor) * math.cos(2 * PI * (day + 10) / 365)
        )
        self.sine_solar_height_season_offset = math.sin(lat_rad) * math.sin(
            self.sun_declination
        )
        self.sine_solar_height_amplitude = math.cos(lat_rad) * math.cos(
            self.sun_declination
        )
        # qili: aob is intermediate variable
        aob = self.sine_solar_height_season_offset / self.sine_solar_height_amplitude
        # qili: astronomical daylength
        self.astro_daylength = 12 * (1 + 2 * math.asin(aob) / PI)
        self.total_sine_solar_elevation = 3600 * (
            self.astro_daylength * self.sine_solar_height_season_offset
            + 24 * self.sine_solar_height_amplitude * math.sqrt(1 - aob ** 2) / PI
        )
        self.total_effective_solar_elevation = 3600 * (
            self.astro_daylength
            * (
                self.sine_solar_height_season_offset
                + 0.4
                * (
                    self.sine_solar_height_season_offset ** 2
                    + 0.5 * self.sine_solar_height_amplitude ** 2
                )
            )
            + 12
            * self.sine_solar_height_amplitude
            * (2.0 + 3.0 * 0.4 * self.sine_solar_height_season_offset)
            * math.sqrt(1 - aob ** 2)
            / PI
        )
        self.sol_const = 1370 * (1 + 0.033 * math.cos(2 * PI * day / 365))
        # qili: tiny diff to original fortran code
        print("solar constant: ", self.sol_const) if self.dbg else None

    def greenhouse_trans_interp(self):
        """
        Interpolation procedure for calculation of transmission of greenhouse

        Output:
        self.dir_construct_trans: transmission of the construction for direct radiation
        self.dir_glass_trans: transmission of the glass for direct radiation
        """
        RAD = 0.017453292
        # adaptation of azimuth depending on orientation of greenhouse
        self.solar_azimuth = self.solar_azimuth - self.deviation_orient * RAD

        # conversion to degrees
        A1 = self.solar_azimuth / RAD
        A1 = math.fmod(A1, 180)
        E = self.solar_elevation / RAD

        # if necessary, mirroring of azimuth of sun
        if 180 >= A1 >= 90:
            A = 180 - A1
        if 0 > A1 > -90:
            A = -A1
        if -90 >= A1 >= -180:
            A = 180 + A1
        if 90 > A1 >= 0:
            A = A1

        for i in range(self.num_layer):
            if E < self.elevation_layer[i]:
                ix_min = max(i, 1)
                ix_max = i + 1
                break
        if i == self.num_layer - 1:
            ix_min = self.num_layer
            ix_max = self.num_layer

        # interpolation in azimuth
        construction_trans_1 = onev(
            A,
            self.construction_transmission[ix_min - 1],
            self.azimuth[ix_min - 1],
            self.num_entries[ix_min - 1],
        )
        glass_trans_1 = onev(
            A,
            self.glass_transmission[ix_min - 1],
            self.azimuth[ix_min - 1],
            self.num_entries[ix_min - 1],
        )
        if ix_min == ix_max:
            self.dir_construct_trans = construction_trans_1
            self.dir_glass_trans = glass_trans_1
        else:
            construction_trans_2 = onev(
                A,
                self.construction_transmission[ix_max - 1],
                self.azimuth[ix_max - 1],
                self.num_entries[ix_max - 1],
            )
            glass_trans_2 = onev(
                A,
                self.glass_transmission[ix_max - 1],
                self.azimuth[ix_max - 1],
                self.num_entries[ix_max - 1],
            )
            self.dir_construct_trans = construction_trans_1 + (
                construction_trans_2 - construction_trans_1
            ) * (E - self.elevation_layer[ix_min - 1]) / (
                self.elevation_layer[ix_max - 1] - self.elevation_layer[ix_min - 1]
            )

            self.dir_glass_trans = glass_trans_1 + (glass_trans_2 - glass_trans_1) * (
                E - self.elevation_layer[ix_min - 1]
            ) / (self.elevation_layer[ix_max - 1] - self.elevation_layer[ix_min - 1])

    def calc_sun_position(self):
        """
        Calculation of position of sun at given day of year, time of day and latitude
        """
        RD = PI / 180
        self.sine_solar_elevation = (
            self.sine_solar_height_season_offset
            + self.sine_solar_height_amplitude
            * math.cos(2 * PI * (self.solar_time + 12) / 24)
        )
        self.solar_elevation = math.asin(self.sine_solar_elevation)

        cos_azimuth = -(
            math.sin(self.sun_declination * RD)
            - math.sin(LAT * RD) * self.sine_solar_elevation
        ) / (math.cos(LAT * RD) * math.cos(self.solar_elevation))

        if cos_azimuth < -1:
            cos_azimuth = -1
        elif cos_azimuth > 1:
            cos_azimuth = 1

        self.solar_azimuth = math.acos(cos_azimuth)
        if self.solar_time <= 12:
            self.solar_azimuth = -1 * self.solar_azimuth

        if self.sine_solar_elevation < 0:
            self.sine_solar_elevation = 0

        return

    def calc_light_response_curve(self, leaf_temp, air_co2):
        """
        Determination of EFF-  and PGMAX- values in negative exponential Determination of EFF and PGMAX values in negative exponential.
        Descriptive formulae are used to calculated initial slope and light-saturation value of negative-exponential light response curve. Formulae are partly developed based on theory of Farquhar, von Caemmerer and Berry (1980).

        Parameters:
        EFF0: Potential light use efficiency in absence of oxygen (mg CO2 J-1).
        RM: the mesophyll resistance to CO2 assimilation.
        RS: stomatal resistance to H2O diffusion.
        RB: boundary layer resistance to H20 diffusion.
        RD20: dark respiration at 20 oC.
        GMT: Table for temperature dependence of mesophyll conductance.
        PMMT: Table for temperature dependence of maximal endogenous.
        """
        EFF0 = 0.017
        RS = 50
        RB = 100
        RD20 = 0.05
        Q10RD = 2.0
        # Table for temperature dependence of mesophyll conductance.
        GMT = [0, 0, 5, 0, 15, 0.004, 25, 0.004, 40, 0, 100, 0]

        # Table for temperature dependence of maximal endogenous photosynthetic capacity.
        PMMT = [0, 0, 5, 0, 15, 2, 25, 2, 40, 0, 100, 0]

        # conductance GM is a function of temperature
        GM = linear_interpolation(GMT, 12, leaf_temp)

        # RM is the mesophyll resistance to CO2 assimilation
        if GM < 0.00001:
            RM = 3.0e30
        else:
            RM = 1 / GM

        # Endogenous photosynthetic capacity PMM (mg CO2 m-2 s-1) is a function of temperaure
        PMM = linear_interpolation(PMMT, 12, leaf_temp)

        # CO2 compensation point increases with temperature according to Brooks & Farquhar, 1985
        gamma = 42.7 + 1.68 * (leaf_temp - 25) + 0.012 * (leaf_temp - 25) ** 2

        co2 = max(air_co2, gamma)
        self.init_light_efficiency = EFF0 * (co2 - gamma) / (co2 + 2 * gamma)

        # PNC is maximum as determined by CO2 diffusion, 1.830 is mg CO2 per M3 per vpm, stomatal resistance and boundary layer resistance to CO2 are 1.6 and 1.36 times larger than to water vapour, respectively.
        PNC = (co2 - gamma) * 1.830 / (RM + 1.36 * RB + 1.6 * RS)

        # PNMAX shows saturation with FNC
        if PMM < 0.00001:
            PNMAX = 0
        else:
            PNMAX = min(PNC, PMM)

        # Dark respiration (mg CO2 m-2 s-1)
        RD = RD20 * Q10RD ** (0.1 * (leaf_temp - 20))
        self.light_assimilation_rate = PNMAX + RD

        return

    def calc_assimilation(self):
        """
        This method performs a Gaussian integration over depth of canopy by selecting three different LAI's and computing assimilation at these LAI levels. The integrated variable is FGROS.
        """

        # Gauss weights for three point Gauss
        i_gauss = 3
        x_gauss = [0.1127, 0.5000, 0.8873]
        w_gauss = [0.2778, 0.4444, 0.2778]

        # Prevent math overflow; name change to prevent change of variable value
        sinel = max(0.02, self.sine_solar_elevation)

        # reflection of horizontal and spherical leaf angle distribution
        sqv = math.sqrt(1 - self._scat_coef_leaves)
        horizon_reflect = (1 - sqv) / (1 + sqv)
        sph_reflect = horizon_reflect * 2 / (1 + 1.6 * sinel)

        # extinction coefficient for direct radiation and total direct flux
        kdirbl = (0.5 / sinel) * self._exti_coef_diff_light / (0.8 * sqv)
        kdirt = kdirbl * sqv

        # Section calculating effect of ground refectance of radiation
        # transmissivity T, effective transmissivity TE and effective reflectivity RE for incoming diffuse (1), incoming direct and its diffused components together (2) andc reflected diffuse radiation from the ground surface (3) reckoned in upward direction.
        t1 = math.exp(-self._exti_coef_diff_light * self.LAI_function)
        t2 = math.exp(-kdirt * self.LAI_function)
        t3 = t1
        corr1 = (
            (horizon_reflect - self.ground_ref_coef)
            / (self.ground_ref_coef - 1 / horizon_reflect)
            * t1
            * t1
        )
        corr2 = -sph_reflect * sph_reflect * t2 * t2
        corr3 = -horizon_reflect * horizon_reflect * t3 * t3
        re1 = (horizon_reflect + corr1 / horizon_reflect) / (1 + corr1)
        re2 = (sph_reflect + corr2 / sph_reflect) / (1 + corr2)
        re3 = (horizon_reflect + corr3 / horizon_reflect) / (1 + corr3)
        te1 = (
            t1
            * (horizon_reflect * horizon_reflect - 1)
            / (horizon_reflect * self.ground_ref_coef - 1)
            / (1 + corr1)
        )
        te2 = t2 * (1 - sph_reflect * sph_reflect) / (1 + corr2)
        te3 = t3 * (1 - horizon_reflect * horizon_reflect) / (1 + corr3)

        # Reflected diffused flux at ground surface originating from direct radiation, including secondary reflection.
        phiu = (
            self.ground_ref_coef
            * self.instant_dir_rad_flux
            * te2
            / (1 - re3 * self.ground_ref_coef)
        )

        # selectiohn of depth of canopy, canopy assimilation is set to zero
        self.instant_canopy_assimilation = 0
        for i in range(i_gauss):
            c_LAI = self.LAI_function * x_gauss[i]

            # absorbed fluxed per unit leaf area: diffuse flux, total direct flux, direct component of direct flux.

            visdf = (
                (1 - horizon_reflect)
                * self._exti_coef_diff_light
                * (
                    self.instant_diff_rad_flux
                    * (
                        math.exp(-self._exti_coef_diff_light * c_LAI)
                        + corr1
                        * math.exp(self._exti_coef_diff_light * c_LAI)
                        / horizon_reflect
                    )
                    / (1 + corr1)
                    + phiu
                    * (
                        math.exp(
                            self._exti_coef_diff_light * (c_LAI - self.LAI_function)
                        )
                        + corr3
                        * math.exp(
                            self._exti_coef_diff_light * (self.LAI_function - c_LAI)
                        )
                        / horizon_reflect
                    )
                    / (1 + corr3)
                )
            )

            vist = (
                (1 - sph_reflect)
                * self.instant_dir_rad_flux
                * kdirt
                * (
                    math.exp(-kdirt * c_LAI)
                    + corr2 * math.exp(kdirt * c_LAI) / sph_reflect
                )
                / (1 + corr2)
            )

            visd = (
                (1 - self._scat_coef_leaves)
                * self.instant_dir_rad_flux
                * kdirbl
                * math.exp(-kdirbl * c_LAI)
            )

            # absorbed flux (J/M2 leaf/s) for shaded leaves and assimilation of shaded leaves.
            visshd = visdf + (vist - visd)
            fgrsh = self.light_assimilation_rate * (
                1
                - math.exp(
                    -visshd * self.init_light_efficiency / self.light_assimilation_rate
                )
            )

            # direct flux absorbed by leaves perpendicular on direct beam and assimilation of sunlit leaf area
            vispp = (1 - self._scat_coef_leaves) * self.instant_dir_rad_flux / sinel
            fgrsun = 0
            for j in range(i_gauss):
                vissun = visshd + vispp * x_gauss[j]
                fgrs = self.light_assimilation_rate * (
                    1
                    - math.exp(
                        -vissun
                        * self.init_light_efficiency
                        / self.light_assimilation_rate
                    )
                )
                fgrsun = fgrsun + fgrs * w_gauss[j]

            # fraction sunlit leaf area (frac_sunlit_leaf_area) and local assimilation rate (fgl)
            frac_sunlit_leaf_area = math.exp(-kdirbl * c_LAI)
            fgl = frac_sunlit_leaf_area * fgrsun + (1 - frac_sunlit_leaf_area) * fgrsh

            # integration of local assimilation rate to canopy assimilation (FGROS)
            self.instant_canopy_assimilation += fgl * w_gauss[i]

        self.instant_canopy_assimilation *= self.LAI_function
        return

    def total_gross_assimilation(self):
        # parameter for conversion for radians to degrees (=PI/180)
        RAD = 0.017453292
        self.sun_declination /= RAD
        LONG = 5.7  # qili: unknown magic number

        # assimilation set to zero and three different times of the day (HOUR)
        self.daily_greenhouse_diffuse_rad = 0
        self.daily_greenhouse_direct_rad = 0
        self.daily_gross_assimilation = 0

        SGSOM = 0  # qili: unknown but an intermediate variable

        for i in range(1, 25):
            self.solar_time = i - (1 - LONG / 15) + solar_time_deviation(self.day)
            self.calc_sun_position()
            if self.sine_solar_elevation >= 0.001:
                self.greenhouse_trans_interp()

                # qili: unknow intermediate variable
                # calculate global radiation outside qili: why i-2
                SG = self._current_day_radiation[i - 2]
                SGSOM = SGSOM + SG * 3600
                frac_diff_rad = calc_fraction_diffuse(
                    self.sol_const, SG, self.sine_solar_elevation
                )
                radiation = 0.47 * SG
                diffuse_radiation = frac_diff_rad * radiation
                direct_radiation = (1 - frac_diff_rad) * radiation

                # light climate inside the greenhouse, multiplication with greenhouse
                self.instant_diff_rad_flux = (
                    diffuse_radiation
                    * self.diffuse_rad_transmission
                    * self.cor_transmissivity
                )
                self.direct_rad_transmission = (
                    self.dir_construct_trans * self.dir_glass_trans
                )
                self.instant_dir_rad_flux = (
                    direct_radiation
                    * self.direct_rad_transmission
                    * self.cor_transmissivity
                )
                temperature = self._current_day_temperature[i - 1]
                co2 = self._current_day_co2[i - 1]
                self.calc_light_response_curve(temperature, co2)
                self.calc_assimilation()
            else:
                self.instant_canopy_assimilation = 0
                self.instant_diff_rad_flux = 0
                self.instant_dir_rad_flux = 0

            # integration of assimilation rate to a daily total
            self.daily_gross_assimilation += self.instant_canopy_assimilation * 3.6

            # integration radiation inside greenhouse to a daily total
            self.daily_greenhouse_diffuse_rad += self.instant_diff_rad_flux * 3600
            self.daily_greenhouse_direct_rad += self.instant_dir_rad_flux * 3600

        return

    def truss_growth(self, writer):
        """
        This method calculates SINK STRENGTH of VEGETATIVE UNITS and individual FRUIT TRUSSES in tomato
        """
        # Coefficients for Richards function (fit of weight/growth period as a function of development stage) for 17, 21 and 25 oC averaged Phytotron experiment nr. 18. (Heuvelink & Marcelis, 1989).
        A = 0.1541 / 1.12 * 1.03
        B = 4.3435
        C = 0.2782
        D = 1.3065

        # Factor to bring our values to those of A.N.M. de Koning
        sink_factor = 1.43

        # Growth rates of individual trusses are calculated from anthesis until harvest. STTRU = start of truss growth FR = flowering rate (trusses/day).
        if self.day < self._start_day_of_truss_growth:
            flowering_rate = 0
        else:
            flowering_rate = (
                -0.2863 + 0.1454 * math.log(self.corrected_temperature)
            ) * self.cor_flower_rate

        self._initial_truss += flowering_rate
        self.truss = int(self._initial_truss)

        # Calculation of sinkstrength of individual trusses if self.IAGE(I) < 1 truss was at anthesis today -> sinkstrength = 0 if ITG(I) = 1 truss was already harvested -> sinkstrength = 0

        self.total = 0
        self.total_veg = 0

        if (self.truss < 1) or ((self.truss == 1) and (self.IAGE[0] < 1)):
            self.frac_dry_matter_leaves = self._frac_dry_weight_leaves
            self.frac_dry_matter_stems = 1 - self._frac_dry_weight_leaves
            self.frac_dry_matter_organs = 0
            if self.truss == 1:
                self.IAGE[0] = self.IAGE[0] + 1
        else:
            self.num_veg_units = self.truss + 3
            for i in range(self.truss):
                if (self.if_truss_harvested[i] == 1) or (self.IAGE[i] < 1):
                    self.sink[i] = 0
                else:
                    self.dev_stages_per_truss[i] = (
                        self.dev_stages_per_truss[i]
                        + (
                            1.814
                            + math.log(self.corrected_temperature / 20)
                            * (
                                3.92
                                - 0.2127 * self.dev_stages_per_truss[i]
                                + 0.004505 * self.dev_stages_per_truss[i] ** 2
                                - 0.000024 * self.dev_stages_per_truss[i] ** 3
                            )
                        )
                        * self.cor_fruit_dev
                    )
                    self.temp_array[i] += self.corrected_temperature
                    self.avtemp_truss[i] = self.temp_array[i] / self.IAGE[i]

                    # Calculation of truss sink strength
                    Y = A * (
                        1 + math.exp(-B * (self.dev_stages_per_truss[i] / 100 - C))
                    ) ** (1 / (1 - D))
                    if (i + 1) >= self.truss_pos:
                        self.sink[i] = (
                            self.new_plant_density
                            * self.fruit_per_truss[i]
                            * sink_factor
                            * self._sink_factor
                            * Y
                            * B
                            / (D - 1)
                            / (
                                math.exp(B * (self.dev_stages_per_truss[i] / 100 - C))
                                + 1
                            )
                        )
                    else:
                        self.sink[i] = (
                            self._plant_density
                            * self.fruits_per_truss_table[i]
                            * sink_factor
                            * self._sink_factor
                            * Y
                            * B
                            / (D - 1)
                            / (
                                math.exp(B * (self.dev_stages_per_truss[i] / 100 - C))
                                + 1
                            )
                        )

                self.total += self.sink[i]
                # print("total before:", self.total)

            # print("self.IAGE:", self.IAGE[0])
            # time.sleep(5)
            if self.IAGE[0] == 1:
                self.veg_dev[0] = 40
                self.veg_dev[1] = 24
                self.veg_dev[2] = 12

            for i in range(self.num_veg_units):
                # print("veg_dev[i]:", veg_dev[i])
                # time.sleep(5)
                if (self.itg_veg[i] == 1) or (self.veg_dev[i] >= 100):
                    self.veg_sink[i] = 0
                else:
                    self.veg_dev[i] = (
                        self.veg_dev[i]
                        + 1.814
                        + math.log(self.corrected_temperature / 20.0)
                        * (
                            3.92
                            - 0.2127 * self.veg_dev[i]
                            + 0.004505 * self.veg_dev[i] ** 2
                            - 0.000024 * self.veg_dev[i] ** 3
                        )
                    )
                    # Calculation of sink strength of vegetative units
                    y_veg = A * (1 + math.exp(-B * (self.veg_dev[i] / 100 - C))) ** (
                        1 / (1 - D)
                    )
                    # print("y_veg: ", veg_dev[i])
                    # time.sleep(5)
                    if (i + 1) >= self.veg_pos:
                        self.veg_sink[i] = (
                            self.new_plant_density
                            * 3
                            * self.cor_veg_sink
                            * sink_factor
                            * self._sink_factor
                            * y_veg
                            * B
                            / (D - 1)
                            / (math.exp(B * (self.veg_dev[i] / 100 - C)) + 1)
                        )
                    else:
                        self.veg_sink[i] = (
                            self._plant_density
                            * 3
                            * self.cor_veg_sink
                            * sink_factor
                            * self._sink_factor
                            * y_veg
                            * B
                            / (D - 1)
                            / (math.exp(B * (self.veg_dev[i] / 100 - C)) + 1)
                        )

                    if i == 0:
                        self.veg_sink[i] *= 2.5

                # TOTAL = total sink strength of VEGETATIVE plant part
                self.total_veg += self.veg_sink[i]

            # Calculation of vegetative sink strength and total plant sink strength
            self.total += self.total_veg
            # print("total:", self.total)
            # time.sleep(10)

            for i in range(self.truss):
                # IF dev_stages_per_truss(I).GE.100. if_truss_harvested(I) = 1; truss removed from plant Write number of truss (I), average experienced temperature and fruit growth period to terminal and to file FruitDevInfo

                if self._manual_leaf_picking:
                    if i == 1:
                        self.num_picked_units = int(
                            linear_interpolation(
                                self._PDVU, self.PDVU_table_size, self.day
                            )
                        )
                        for j in range(self.num_picked_units):
                            self.itg_veg[j] = 1
                else:
                    if self.dev_stages_per_truss[i] >= (90 * self.cor_leaf_picking):
                        self.itg_veg[i] = 1

                if (self.dev_stages_per_truss[i] >= 100) and (self.if_truss_harvested[i] != 1):
                    self.if_truss_harvested[i] = 1

                    if self.dbg:
                        print(
                            f"Truss number:{i}, day: {self.day}, avtemp_truss: {self.avtemp_truss[i]}, number of days from anthesis to harvest: {self.IAGE[i]}"
                        )
                    writer.writerow([i + 1, self.day, self.avtemp_truss[i], self.IAGE[i]])
                else:
                    self.IAGE[i] += 1

            self.frac_dry_matter_leaves = (
                self.total_veg / self.total * self._frac_dry_weight_leaves
            )
            self.frac_dry_matter_stems = (
                self.total_veg / self.total * (1 - self._frac_dry_weight_leaves)
            )
            self.frac_dry_matter_organs = (self.total - self.total_veg) / self.total

        return

    def crop_growth(self):
        i_time = int(self.time)
        i_start_time = int(self._start_time)

        # Maintenance respiration [g CH2O m-2 day-1]
        maints = self.cor_maintenance_respiration * (
            self.dry_weight_leaves_existing * self._maint_respiration_coef_leaves
            + self.dry_weight_stems_existing * self._maint_respiration_coef_stems
            + self.dry_weight_organs_existing * self._maint_respiration_coef_organs
            + self._init_dry_weight_roots * self._maint_respiration_coef_roots
        )
        teff = self._Q10_value ** ((self.corrected_temperature - self._ref_temp) / 10)

        if i_time == i_start_time:
            self.maint = min(self.daily_total_photosyn, maints * teff)
            # Assimilate requirements for dry matter conversion [g CH2O/g dry matter]
            assimilate_req = self.cor_growth_respiration * (
                self.frac_dry_matter_leaves * self._assimilate_req_leaves
                + self.frac_dry_matter_stems * self._assimilate_req_stem
                + self.frac_dry_matter_organs * self._assimilate_req_organs
                + self.frac_dry_matter_roots * self._assimilate_req_roots
            )
            self.total_growth_rate = (
                self.daily_total_photosyn - self.maint
            ) / assimilate_req

        self.growth_rate_table[i_time - 1] = self.total_growth_rate / (
            self.dry_weight_leaves_existing
            + self._init_dry_weight_stem
            + self.dry_weight_organs_existing
            + self._init_dry_weight_roots
        )
        if i_time >= (i_start_time + 4):
            self.av_growth_rate = (
                self.growth_rate_table[i_time - 1]
                + self.growth_rate_table[i_time - 2]
                + self.growth_rate_table[i_time - 3]
                + self.growth_rate_table[i_time - 4]
                + self.growth_rate_table[i_time - 5]
            ) / 5
        elif i_time == i_start_time:
            self.av_growth_rate = self.growth_rate_table[i_time - 1]
        else:
            self.RGRSOM += self.growth_rate_table[i_time - 1]
            self.av_growth_rate = self.RGRSOM / (int(self.time - self._start_time))

        maint_help = (
            maints * teff * self.A * (1 - math.exp(-self.B * self.av_growth_rate))
        )
        self.maint = min(self.daily_total_photosyn, maint_help)

        # Assimilate requirements for dry matter conversion [g CH2O/g dry matter]
        assimilate_req = self.cor_growth_respiration * (
            self.frac_dry_matter_leaves * self._assimilate_req_leaves
            + self.frac_dry_matter_stems * self._assimilate_req_stem
            + self.frac_dry_matter_organs * self._assimilate_req_organs
            + self.frac_dry_matter_roots * self._assimilate_req_roots
        )
        
        self.total_growth_rate = (
            self.daily_total_photosyn - self.maint
        ) / assimilate_req


        self.total_growth_rate += self.buffer
        self.total = self.total_veg * (1 + self.fruit_veg_ratio) + (
            self.total - self.total_veg
        )

        # Calculation of growth of vegetative plant parts
        self.truss = int(self._initial_truss)
        
        if (self.truss < 1) or (self.truss == 1) and (self.IAGE[0] <= 1):
            self.leaves_growth = self.frac_dry_matter_leaves * self.total_growth_rate
            self.stems_growth = self.frac_dry_matter_stems * self.total_growth_rate
            self.roots_growth = self.frac_dry_matter_roots * self.total_growth_rate
            self.fruits_growth = 0
        else:
            for i in range(self.num_veg_units):
                if self.total >= self.total_growth_rate:
                    self.veg_growth[i] = (
                        self.veg_sink[i] / self.total * self.total_growth_rate
                    )
                else:
                    self.veg_growth[i] = self.veg_sink[i]

            if self.total >= self.total_growth_rate:
                self.roots_growth = self.total_growth_rate * self.frac_dry_matter_roots
                self.leaves_growth = (
                    self.total_growth_rate * self.frac_dry_matter_leaves
                )
                self.stems_growth = self.total_growth_rate * self.frac_dry_matter_stems
                self.buffer = 0.0
            else:
                self.roots_growth = (
                    self.total_veg / (1 - self.fruit_veg_ratio) * self.fruit_veg_ratio
                )
                self.leaves_growth = self.total_veg * self._frac_dry_weight_leaves
                self.stems_growth = self.total_veg * (1 - self._frac_dry_weight_leaves)
                self.buffer = self.total_growth_rate - self.total

            # Calculation of growth of individual trusses
            self.fruits_growth = 0
            for i in range(self.truss):
                if self.total >= self.total_growth_rate:
                    self.growth_truss[i] = (
                        self.sink[i] / self.total * self.total_growth_rate
                    )
                else:
                    self.growth_truss[i] = self.sink[i]
                self.fruits_growth += self.growth_truss[i]

        if not self._fixed_crop:
            for i in range(self.truss):
                self.truss_weight_table[i] += self.growth_truss[i]
            for i in range(self.num_veg_units):
                self.veg_dry_weight_unit[i] += self.veg_growth[i]

            self._init_dry_weight_roots += self.roots_growth
            self._init_dry_weight_leaves += self.leaves_growth
            self._init_dry_weight_stem += self.stems_growth
            self._init_dry_weight_orgrans += self.fruits_growth
            self.total_growth_rate = (
                self.leaves_growth
                + self.stems_growth
                + self.fruits_growth
                + self.roots_growth
            )
            self.dry_weight_stems_existing = self._init_dry_weight_stem
            self.dry_weight_organs_existing = 0

            for i in range(self.truss):
                if self.if_truss_harvested[i] != 1:
                    self.dry_weight_organs_existing += self.truss_weight_table[i]

            weight_leaves_picked = 0
            if not self._manual_leaf_picking:
                if (self.dev_stages_per_truss[0] >= (70 * self.cor_leaf_picking)) and (
                    self.dev_stages_per_truss[0] < (90 * self.cor_leaf_picking)
                ):
                    weight_leaves_picked = self.veg_dry_weight_unit[0] * 0.7 / 2
                else:
                    for i in range(self.num_veg_units):
                        if self.itg_veg[i] == 1:
                            weight_leaves_picked += self.veg_dry_weight_unit[i] * 0.7
            else:
                for i in range(self.num_veg_units):
                    if self.itg_veg[i] == 1:
                        weight_leaves_picked += self.veg_dry_weight_unit[i] * 0.7

            self.num_picked_units = 0
            for i in range(self.num_veg_units):
                if self.itg_veg[i] == 1:
                    self.num_picked_units += 1

            self.dry_weight_leaves_existing = (
                self._init_dry_weight_leaves - weight_leaves_picked
            )

            # qili: a useless statement here
            # HARVEST = HARVEST / 55
            self.total_dry_weight = (
                self._init_dry_weight_leaves
                + self._init_dry_weight_stem
                + self._init_dry_weight_orgrans
                + self._init_dry_weight_roots
            )
        return

    def print_growth(self, growth_writer, fruit_dev_writer):

        if self._start_time == self.time:
            growth_header = [
                "time (jday)",
                "daily_total_rad (MJ/m2/d)",
                "temperature (Celsius)",
                "CO2",
                "light_assimilation_rate (mgCO2/m2/s)",
                "init_light_efficiency (mgCO2/J)",
                "SLA*10000",
                "LAI",
                "daily total gross photosynthesis (g CH20/m2/d)",
                "average_growth_rate (1/d)",
                "maintenance respiration (gCH2O/m2/d)",
                "total_growth_rate (g/m2/d)",
                "dry_weight (total) (g/m2)",
                "dry_weight (total - root) (g/m2)",
                "dry_weight (leave) (g/m2)",
                "dry_weight (leave + stem) (g/m2)",
                "dry_weight (organ) (g/m2)",
                "dry_weight (root) (g/m2)",
                "excluding_dry_weight (leave) (g/m2)",
                "excluding_dry_weight (stem) (g/m2)",
                "excluding_dry_weight (organ) (g/m2)",
                "daily_greenhouse_diffuse_radiation (MJ/m2/d)",
                "daily_greenhouse_direct_radiation (MJ/m2/d)",
            ]
            growth_writer.writerow(growth_header)

            fruit_dev_header = [
                "truss_number",
                "day",
                "average_temperature",
                " number of days from anthesis to harvest",
            ]
            fruit_dev_writer.writerow(fruit_dev_header)

        if (
            (self.time == self._start_time)
            or (self.icount == self._output_time_step)
            or (self.time == self._finish_time)
        ):
            data_row = [
                self.time,
                self.daily_total_rad / 1000000,
                self.corrected_temperature,
                self.corrected_co2,
                self.light_assimilation_rate,
                self.init_light_efficiency,
                self.SLA_function * 10000,
                self.LAI_function,
                self.daily_total_photosyn,
                self.av_growth_rate,
                self.maint,
                self.total_growth_rate,
                self.total_dry_weight,
                self.total_dry_weight - self._init_dry_weight_roots,
                self._init_dry_weight_leaves,
                self._init_dry_weight_leaves + self._init_dry_weight_stem,
                self._init_dry_weight_orgrans,
                self._init_dry_weight_roots,
                self.dry_weight_leaves_existing,
                self.dry_weight_stems_existing,
                self.dry_weight_organs_existing,
                self.daily_greenhouse_diffuse_rad / 1e6,
                self.daily_greenhouse_direct_rad / 1e6,
            ]
            for i in range(len(data_row)):
                data_row[i] = round(data_row[i], 3)
            self.icount = 0
            growth_writer.writerow(data_row)

        self.icount += 1
        return

    def print_distribution(self, writer):
        if self._start_time == self.time:
            row_header = [
                "time (jday)",
                "fruit_dry_matter_fraction",
                "total_veg_dry_matter",
                "total_dry_matter",
                "SLA",
                "growth_rate",
                "dry_weight (total)",
                "dry_weight (leaves)",
                "dry_weight (stems)",
                "dry_weight (organs)",
                "buffer",
                "truss_number",
                "temperature",
                "CO2",
                "veg_unit_number",
            ]
            for i in range(1, 36):
                row_header.append("TRUSS_" + str(i))
            row_header.append("WTRUSS_615")
            for i in range(1, 17):
                row_header.append("WVEGUNIT_" + str(i))
            writer.writerow(row_header)

        if (
            (self.time == self._start_time)
            or (self.icount1 == self._output_time_step)
            or (self.time == self._finish_time)
        ):
            data_row = [
                self.time,
                (self.frac_dry_matter_leaves + self.frac_dry_matter_stems)
                / (
                    self.frac_dry_matter_leaves
                    + self.frac_dry_matter_stems
                    + self.frac_dry_matter_organs
                ),
                self.total_veg,
                self.total,
                self.SLA_function * 10000,
                self.growth_rate,
                self.total_dry_weight,
                self._init_dry_weight_leaves,
                self._init_dry_weight_stem,
                self._init_dry_weight_orgrans,
                self.buffer,
                self._initial_truss,
                self.corrected_temperature,
                self.corrected_co2,
                self.num_veg_units,
            ]
            for i in range(35):
                data_row.append(self.truss_weight_table[i])
            data_row.append(self.wtruss_615)
            for i in range(16):
                data_row.append(self.veg_dry_weight_unit[i])
            self.icount1 = 0
            # decimal accuracy
            for i in range(len(data_row)):
                data_row[i] = round(data_row[i], 3)
            writer.writerow(data_row)

        self.icount1 += 1
        return

    def start_simulation(self, config_path):
        self.read_radiation()
        self.read_config(config_path)
        self.read_greenhouse_config(config_path)
        self.pre_simulation()

        print("simulation start...") if self.dbg else None
        output_growth = "output/tomsim_growth.csv"
        output_distri = "output/tomsim_distri.csv"
        output_FruitDev = "output/tomsim_FruitDev.csv"

        f_growth = open(output_growth, "w", encoding="utf-8", newline="")
        f_distri = open(output_distri, "w", encoding="utf-8", newline="")
        f_fruitdev = open(output_FruitDev, "w", encoding="utf-8", newline="")

        growth_writer = csv.writer(f_growth)
        distri_writer = csv.writer(f_distri)
        fruitdev_writer = csv.writer(f_fruitdev)
        self.truss_counter = []
        while self.time <= self._finish_time:
            self.itime = self.time

            # calculate Julian day number from self.time
            self.day = (self.time - 1) % 365 + 1
            self.daily_total_rad = 0
            self.cor_daily_total_rad = (
                self._daily_radiation_sum_table[self.itime - 1] * self.cor_irradiance
            )
            for i in range(24):
                self._current_day_radiation[i] = (
                    self._hourly_radiation_table[self.itime - 1, i]
                    * self.cor_irradiance
                )
                self._current_day_temperature[i] = (
                    self._hourly_temperature_table[self.itime - 1, i]
                    * self.cor_temperature
                )
                self._current_day_co2[i] = (
                    self._hourly_co2_table[self.itime - 1, i] * self.cor_co2
                )
                self.daily_total_rad += self._current_day_radiation[i] * 3600

            self.corrected_temperature = (
                self._daily_average_temperature_table[self.itime - 1]
                * self.cor_temperature
            )
            self.corrected_co2 = self._daily_average_co2_table[self.itime - 1]
            print(
                "temp and co2:\n", self.corrected_temperature, self.corrected_co2
            ) if self.dbg else None

            if self._calc_LAI_from_SLA:
                if self._calc_SLA_standard:
                    self.SLA_function = (
                        266 + 88 * math.sin(2 * PI * (self.day + 68) / 365)
                    ) / 10000
                else:
                    self.SLA_function = (
                        linear_interpolation(
                            self._SLA_table, self.SLA_table_size, self.time
                        )
                        / 10000
                    )
                self.SLA_function *= self.cor_SLA
                self.LAI_function = self.SLA_function * self.dry_weight_leaves_existing
            else:
                self.LAI_function = linear_interpolation(
                    self._LAI_table, self.LAI_table_size, self.time
                )

            # calculate dry matter distribution ratios
            self.fruit_veg_ratio = 0.15
            self.lat_radiation = LAT * PI / 180
            self.daytime_radiation_character(self.day, self.lat_radiation)
            self.total_gross_assimilation()

            # Gross photosynthesis [g CH2O m-2 day-1]
            self.daily_total_photosyn = self.daily_gross_assimilation * 30 / 44

            if self._side_shoot:
                if (self.day >= self._flower_shoot_date) and (not self._initiation):
                    self.truss_pos = self._initial_truss
                    self.veg_pos = self.num_veg_units - 3
                    self._initiation = True

            # PT2000 CORFR/CORDEV added to call for sensitivity analysis
            self.truss_growth(fruitdev_writer)

            fractions = (self.frac_dry_matter_leaves + self.frac_dry_matter_stems) * (
                1 + self.fruit_veg_ratio
            ) + self.frac_dry_matter_organs

            self.frac_dry_matter_roots = (
                self.fruit_veg_ratio
                * (self.frac_dry_matter_leaves + self.frac_dry_matter_stems)
                / fractions
            )

            self.frac_dry_matter_leaves /= fractions
            self.frac_dry_matter_stems /= fractions
            self.frac_dry_matter_organs /= fractions

            self.crop_growth()

            for i in range(5, int(self._initial_truss)):
                self.wtruss_615 += self.growth_truss[i]

            self.print_growth(growth_writer, fruitdev_writer)

            self.growth_rate = self.total_growth_rate
            self.print_distribution(distri_writer)
            self.truss_counter.append(self._initial_truss)
            self.time += self._sim_time_step

        f_growth.close()
        f_distri.close()

        community_1 = 0
        for i in range(30):
            community_1 += self.truss_weight_table[i] / 30
        community_1 = (
            community_1
            / self._plant_density
            / self.fruit_per_truss
            / self.cor_fruit_per_truss
        )
        fruitdev_writer.writerow(
            ["average final truss dry weight 1-30 ", round(community_1, 3)]
        )
        i_truss = int(self._initial_truss)
        community_2 = 0
        fruitdev_writer.writerow(["final truss number: ", i_truss])

        for i in range(i_truss):
            community_2 += self.truss_weight_table[i] / i_truss
        community_2 = (
            community_2
            / self._plant_density
            / self.fruit_per_truss
            / self.cor_fruit_per_truss
        )
        fruitdev_writer.writerow(
            ["average truss DW based on truss-number:", community_2]
        )

        f_fruitdev.close()

        total_dry_yield = self.total_dry_weight - (self._init_dry_weight_roots + self.dry_weight_stems_existing + self.dry_weight_leaves_existing + self.dry_weight_organs_existing)

        days = [*range(1, self._finish_time - self._start_time + 2)]
        print(days)
        print(self.truss_counter)
        # the dry weight distribution pie chart data
        dry_weight_distribution = {
            "leaves":self.dry_weight_leaves_existing / self.total_dry_weight,
            "organs":self.dry_weight_organs_existing / self.total_dry_weight,
            "stems":self.dry_weight_stems_existing / self.total_dry_weight,
            "roots":self._init_dry_weight_roots / self.total_dry_weight,
            "yield":total_dry_yield / self.total_dry_weight
        }
        truss_growh = {
            "days":days,
            "truss_number": self.truss_counter
        }
        # the tomatos contain approximately %94 water.
        fresh_yield = total_dry_yield / 0.06

        print("simulation is over...") if self.dbg else None

        return fresh_yield, dry_weight_distribution, truss_growh