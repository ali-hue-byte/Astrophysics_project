# Astrophysics_Project


## Description
This project analyzes astrophysical data and simulates the motion of celestial bodies using numerical methods.

It was developed as part of an Imperative Programming course in a team of two students using Python. The project includes:
- Statistical analysis of real astrophysical datasets
- Verification of physical laws (Kepler’s laws, power-law relations)
- A 2D N-body simulation of the Solar System using the Verlet integration method

## Features

### Part 1: Analysis
- **Data analysis with linear regression**: We used linear regression to find relationships between real astrophysical data and compare them with theoretical predictions. The least squares method was implemented from scratch.
- **Logarithmic transformation**: To study some power-law relationships, we transformed the data to a logarithmic scale so we can apply linear regression.
- **Kepler's third law**: Verified the proportionality between the square of the orbital period and the cube of the semi-major axis using data from Jupiter's moons and Solar System planets.
- **Temperature-distance relationship**: Studied the relationship between exoplanets equilibrium temperature and their distance from their star.
- **Mass-radius relationship of stars**: Analyzed the relationship between the mass and radius of stars using real stellar data, for both stars more and less massive than the Sun.

### Part 2: Simulation
- **Solar System**: Animated simulation of planetary trajectories using Verlet integration method and real NASA ephemeris data for initial positions and velocities.
- **Random new planets**: Added new planets with random characteristics.
- **Asteroids**: Added probability of appearance of randomly generated asteroids.
- **Collision detection**: The program detects collisions between bodies and merges them using the *conservation of linear momentum*.

---

## Technologies

- Python
- NumPy
- Matplotlib
  
---

## Project Structure

### Python files
- **Astrophysics_project_en.py**: Main program (analysis + simulation, English comments)
- **Astrophysics_project_fr.py**: Same program with French comments

### Reports
- **Project_Report_en**: Full explanation of methods, physics laws, and architecture (English)
- **Rapport_Projet_fr**: Same report in French

### Datasets
- **Stars.csv**: Stellar mass and radius data used for mass–radius relation study
- **Solar_System_data.csv**: Orbital data of Solar System planets (Kepler’s law analysis)
- **Jupiter_Moons.csv**: Orbital parameters of Jupiter’s moons
- **Exoplanets.csv**: Exoplanet catalog used for temperature–distance study

---

## Installation

1. Clone the repository:

```bash
git clone https://github.com/ali-hue-byte/Astrophysics_project
cd Astrophysics_project
```
2. Install the required dependencies:

```bash
pip install -r requirements.txt
```
3. Run the main application:
Make sure you run the script from the project root directory so that all datasets can be accessed correctly.

```bash
python Astrophysics_project_en.py
```

or 

```bash
python Astrophysics_project_fr.py
```
---

## Screenshots
### Temperature–Distance Relationship (Exoplanets)
<img width="1559" height="1244" alt="image" src="https://github.com/user-attachments/assets/c2f613b2-866c-4653-be2b-fdb8f76f8133" />

### Solar System N-body Simulation
<img width="2545" height="1376" alt="image" src="https://github.com/user-attachments/assets/83ec4277-c4f4-46bf-9a9e-0ef829678138" />
