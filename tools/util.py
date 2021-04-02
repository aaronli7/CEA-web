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