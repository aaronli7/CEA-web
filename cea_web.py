from flask import Flask, request, render_template, jsonify, redirect
import requests
import sqlite3
import time, os, datetime
from tools.util import *
from CEAsimulator.TomSim import TomSim

app = Flask(__name__, template_folder="templates")
db_path = "us_cities.db"
sql_path = "us_cities.sql"
date_format ="%Y-%m-%d"
simulator_config = os.path.abspath(os.path.join(os.getcwd(), "CEAsimulator/config.ini"))

@app.route("/signup", methods=['POST', 'GET'])
def signup():
    if request.method == 'POST':
        username = request.form["username"]
        nameofUser = request.form["name"]
        password = request.form["password"]
        uploaded_file = request.files['file']
        if uploaded_file.filename != '':
            #filename = secure_filename(uploaded_file.filename)
            uploaded_file.save(uploaded_file.filename)
        with sqlite3.connect(db_path) as con:
            con.execute("INSERT INTO users(username, nameofUser, password) VALUES (?, ?, ?)",(username, nameofUser, password))
            temp = "SELECT * FROM users"
            pdtemp = pd.read_sql(temp, con, index_col=None)
            print(pdtemp)
        return render_template("signupSuccess.html", title="CEA Simulator Signup Successful", name = nameofUser, file = uploaded_file)

    return render_template("signup.html", title="CEA Simulator signup Page")

@app.route("/signupSuccess", methods=['POST', 'GET'])
def signupSuccess():
    if request.method == 'POST':
        
        return redirect("/growthSimulation")
    return render_template("signupSuccess.html", title="CEA Simulator Signup Successful")

@app.route("/growthSimulation", methods=['POST', 'GET'])
def home():
    if request.method == 'POST':
        co2 = float(request.form["CO2"])
        temp = float(request.form["temperature"])
        city_id = int(request.form["city_id"])
        fruit_per_truss = int(request.form["fruit_per_truss"])
        start_date = (request.form["start_date"])
        end_date = (request.form["end_date"])
        start_julian_day = datetime.datetime.strptime(start_date, date_format)
        end_julian_day = datetime.datetime.strptime(end_date, date_format)

        start_julian_day = start_julian_day.timetuple().tm_yday
        end_julian_day = end_julian_day.timetuple().tm_yday
        
        query_location = f"""
            SELECT 
                LATITUDE, LONGITUDE 
            FROM US_CITIES 
            WHERE ID={city_id}
        """
        with sqlite3.connect(db_path) as con:
            location = pd.read_sql(query_location, con, index_col=None)
            lat = location.iloc[0]["LATITUDE"]
            lon = location.iloc[0]["LONGITUDE"]
            # output_array = location.to_json(orient="values")
            sim = TomSim(co2=co2, temperature=temp, fruit_per_truss=fruit_per_truss, lon=lon, lat=lat, start_date=start_julian_day, end_date=end_julian_day, debug=False)
            fresh_yield, dryweight_distribution, truss_growth = sim.start_simulation(config_path=simulator_config)
            
            return jsonify(
                latitude = lat,
                longitude = lon,
                total_yield = fresh_yield,
                leaves = dryweight_distribution['leaves'],
                organs = dryweight_distribution['organs'],
                stems = dryweight_distribution['stems'],
                roots = dryweight_distribution['roots'],                
                yielded = dryweight_distribution['yield']
            )
            

    show_states_query = "SELECT ID, STATE_NAME FROM US_STATES"
    with sqlite3.connect(db_path) as con:
        result = pd.read_sql(show_states_query, con, index_col=None)
        states_name = result['STATE_NAME']
        states_id = result['ID']
        states_info = zip(states_id, states_name)
    return render_template("home.html", title="CEA Simulator Home Page", states_info=states_info)

@app.route("/city", methods=["POST", "GET"])
def select_states():
    if request.method == 'POST':
        data_type = request.form["type"]

        if data_type == "countyData":
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
                
        elif data_type == "cityData":
            city_info = request.form.getlist("category_id[]")
            show_city_query = f"""
                SELECT
                    ID, CITY
                FROM US_CITIES
                WHERE ID_STATE={city_info[0]} AND COUNTY='{city_info[1]}'
                ORDER BY CITY ASC;
            """
            with sqlite3.connect(db_path) as con:
                result = pd.read_sql(show_city_query, con, index_col=None)
                output_array = result.to_json(orient="values")

    return output_array

if not os.path.exists(db_path):
    create_db(db=db_path, sql_script=sql_path)

if __name__ == "__main__":
    app.run(debug=True, port=5000)