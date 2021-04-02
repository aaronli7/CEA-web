import sqlite3
import numpy as np
import pandas as pd
import random
import datetime, time
import sys, requests, io

random.seed(8000)
np.random.seed(8000)

NREL_API_KEY = "qZauaW6FeYkiDGqzYkqfAft9FgFwWJmd3T7UPVEJ"
NREL_EMAIL = "qi.li@uga.edu"
NREL_URL = "https://developer.nrel.gov/api/nsrdb/v2/solar/psm3-tmy-download.csv"

def create_db(db, sql_script):
    """
    Create the database and the tables, including the some test data of 3 resources, users and 4 RFID scanners.
    """
    with sqlite3.connect(db) as con, open(sql_script) as script:
        con.executescript(script.read())
        con.commit()

    return



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

    print(df)
    return df

query_ghi(-161.207778, 55.999722)