import math
import requests
import io
import pandas as pd

NREL_API_KEY = "qZauaW6FeYkiDGqzYkqfAft9FgFwWJmd3T7UPVEJ"
NREL_EMAIL = "qi.li@uga.edu"
NREL_URL = "https://developer.nrel.gov/api/nsrdb/v2/solar/psm3-tmy-download.csv"

def linear_interpolation(table, table_size, x):
    """
    This funciton is a linear interpolation function.
    """
    if table_size % 2 != 0:
        print("ERROR in function: table_size = {}, must be even!".format(table_size))
        quit()
    iup = 0
    for i in range(3, table_size, 2):
        if table[i - 1] <= table[i - 3]:
            print(
                "X-coordinates not in ascending order at element {}, function contains {} points".format(
                    i - 1, table_size
                )
            )
            quit()
        if iup == 0 and table[i - 1] >= x:
            iup = i

    if x < table[0]:
        print(
            "Interpolation below defined region!!\nlint-function contains {} points, Interpolation at X={}".format(
                table_size, x
            )
        )
        return table[1]

    if x > table[table_size - 1 - 1]:
        print(
            "Interpolation above defined region!!\nlint-function contains {} points, Interpolation at X={}".format(
                table_size, x
            )
        )
        return table[table_size - 1]

    # normal interpolation
    slope = (table[iup + 1 - 1] - table[iup - 1 - 1]) / (
        table[iup - 1] - table[iup - 2 - 1]
    )
    return table[iup - 1 - 1] + (x - table[iup - 2 - 1]) * slope


def solar_time_deviation(day):
    """
    Calculation of deviation of solar time as a result of variation of length of solar day during the year because of excentricity of elliptical track of earth around sun, and because of obliquity of ecliptical plane of earth.

    input: day (day of the year)
    output: time deviation (hour)
    """
    PI = 3.1415926
    if 106 >= day >= 1:
        return (-14.2 / 60) * math.sin((day + 7) * PI / 111)
    elif 166 >= day >= 107:
        return (4 / 60) * math.sin((day - 106) * PI / 59)
    elif 246 >= day >= 167:
        return (-6.5 / 60) * math.sin((day - 166) * PI / 80)
    elif 366 >= day >= 247:
        return (16.4 / 60) * math.sin((day - 247) * PI / 113)

    return


def onev(Y, transmission, azimuth, num_of_entry):
    """
    Interpolation in azimuth
    """
    for i in range(num_of_entry):
        if Y < azimuth[i]:
            break
    if i == num_of_entry - 1:
        return transmission[num_of_entry - 1]

    if i == 0:
        return transmission[i]
    else:
        return transmission[i - 1] + (transmission[i] - transmission[i - 1]) * (
            Y - azimuth[i - 1]
        ) / (azimuth[i] - azimuth[i - 1])


def calc_fraction_diffuse(solar_const, global_rad, sine_solar_elevation):
    """
    Calculation of fraction diffuse in global radiation from relation of hourly values of atmosferic transmission versus hourly values of fraction diffuse.

    Input:
    solar_const: solar constant [J m-2 s-1]
    global_rad: global radiation [J m-2 s-1]
    sine_solar_elevation: sine of solar elevation

    Output:
    frac_diff_rad: fraction diffuse radiation
    """

    # qili: unknown intermediate variables
    so = solar_const * sine_solar_elevation
    atmtr = global_rad / so
    frac_diff_rad = 1.47 - 1.66 * atmtr

    if 0.35 >= atmtr > 0.22:
        frac_diff_rad = 1 - 6.4 * (atmtr - 0.22) ** 2
    elif atmtr <= 0.22:
        frac_diff_rad = 1
    frac_diff_rad = max(
        frac_diff_rad, 0.15 + 0.85 * (1 - math.exp(-0.1 / sine_solar_elevation))
    )

    return frac_diff_rad


def query_ghi(lon, lat):
    """
    Query the annual GHI of a given location through NREL API.
    return as a pandas dataframe.
    """
    params = {
        "api_key": NREL_API_KEY,
        "attributes": "ghi",
        "names": "tmy",
        "email": NREL_EMAIL,
        "wkt": f"POINT({lon} {lat})"   
    }

    r = requests.get(url = NREL_URL, params=params)

    data = r.text.split("\n",2)[2] #remove the first two lines
    data = io.StringIO(data)
    df = pd.read_csv(data, sep=",")

    # change to the column Day's format to the what CEAsimulator needs
    day_list = []
    for i in range(1, 366):
        day_list.extend([i] * 24)
    df["Day"] = day_list

    return df