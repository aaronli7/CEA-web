from flask import Flask, request, render_template
import requests
import sqlite3
import time, os
from tools.util import *
from CEAsimulator import TomSim

app = Flask(__name__, template_folder="templates")
db_path = "us_cities.db"
sql_path = "us_cities.sql"


@app.route("/", methods=['POST', 'GET'])
def home():
    if request.method == 'POST':
        pass

    show_states_query = "SELECT ID, STATE_NAME FROM US_STATES"
    with sqlite3.connect(db_path) as con:
        result = pd.read_sql(show_states_query, con, index_col=None)
        states_name = result['STATE_NAME']
        states_id = result['ID']
        states_info = zip(states_id, states_name)
    return render_template("home.html", title="CEA Simulator Home Page", states_info=states_info)


if not os.path.exists(db_path):
    create_db(db=db_path, sql_script=sql_path)

if __name__ == "__main__":
    app.run(debug=True, port=5000)