from flask import Flask, request, render_template, jsonify
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

@app.route("/state", methods=["POST", "GET"])
def select_states():
    if request.method == 'POST':
        state_id = request.form["category_id"]
        show_county_query = f"""
            SELECT
                DISTINCT COUNTY
            FROM US_CITIES
            WHERE ID_STATE={state_id}
            ORDER BY COUNTY ASC;
        """
        with sqlite3.connect(db_path) as con:
            result = pd.read_sql(show_county_query, con, index_col=None)
            county_name = result['COUNTY']
            output_array = county_name.to_json(orient="values")
            # print(output_array[0]["COUNTY"])
    return output_array

# @app.route("/city", methods=["POST", "GET"])
# def select_city():
#     if request.method == 'POST':
#         state_id = request.form["category_id"]
#         show_county_query = f"""
#             SELECT
#                 DISTINCT COUNTY
#             FROM US_CITIES
#             WHERE ID_STATE={state_id}
#         """
#         with sqlite3.connect(db_path) as con:
#             result = pd.read_sql(show_county_query, con, index_col=None)
#             county_name = result['COUNTY']
#             output_array = county_name.to_json(orient="values")
#             # print(output_array[0]["COUNTY"])
#     return output_array

if not os.path.exists(db_path):
    create_db(db=db_path, sql_script=sql_path)

if __name__ == "__main__":
    app.run(debug=True, port=5000)