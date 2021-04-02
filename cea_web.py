from flask import Flask, request, render_template
import requests
import sqlite3
import time, os
from tools.util import *
from CEAsimulator import TomSim

app = Flask(__name__, template_folder="templates")
db_path = "us_cities.db"
sql_path = "us_cities.sql"


@app.route("/")
def home():
    return render_template("home.html", title="CEA Simulator Home Page")


if not os.path.exists(db_path):
    create_db(db=db_path, sql_script=sql_path)

if __name__ == "__main__":
    app.run(debug=True, port=5000)