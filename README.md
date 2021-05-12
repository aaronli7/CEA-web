# INFO 8000 course project (Spring 2021)
## Introduction
**Group menmbers**: *Qi Li*, *Soumya pal*

In this project, we built a web application for users to simulate the tomato growth in the greenhouse.

By inputing the parameters such as temperature, CO2, simulation time range and location, users can get the numerical and visulization result including **total yield of the tomatoes**, **dry weight distribution** and **Growth of truss number**.

Please visit [tomato simulator!](https://info8000.aranisme.top/)

## Components
**database tables**:
- users: storing user information
- cities and states: storing all cities in the United States with coordinates.

**Signup page**: For new user to sign up. User can upload a image. [Signup](https://info8000.aranisme.top/signup)

**Login page**: For user to login. [Login](https://info8000.aranisme.top/login)

**Main page**: For user to simulate the tomato growth and check the result. 

## Tools
**Flask, bootstrap, Jquery, canvasJS, Ajax**

## Deployment
The web application has been deployed on a google cloud platform at: [tomato simulator!](https://info8000.aranisme.top/)

Login:
```
ssh -i ~/.ssh/info8000sp21 kjjohnsen@34.74.130.215
```