# Astrophysics_Project


## Description
Have you ever heard of astrophysics? Probably yes. It is one of the most amazing and challenging fields of science, combining physics and astronomy to study the universe and its celestial bodies.
This project was developed as part of the *Imperative Programming* course at university, in a team of two students. We used Python to analyze real data using linear regression to understand some fundamental laws of astrophysics. We also created a simulation of the Solar System to visualize its behavior and functionality.

## Features

### Part 1: Analysis
- **Data analysis with linear regression**: We used linear regression to find relationships between real astrophysical data and compare them with theoretical predictions. The least squares method was implemented from scratch.
- **Logarithmic transformation**: To study some power-law relationships, we transformed the data to a logarithmic scale so we can apply linear regression.
- **Kepler's third law**: Verified the proportionality between the square of the orbital period and the cube of the semi-major axis using data from Jupiter's moons and solar system planets.
- **Temperature-distance relationship**: Studied the relationship between exoplanets equilibrium temperature and their distance from their star.
- **Mass-radius relationship of stars**: Analyzed the relationship between the mass and radius of stars using real stellar data, for both stars more and less massive than the Sun.

### Part 2: Simulation
- **Solar system simulation**: Animated simulation of the solar system planets trajectories using verlet integration method and real data from NASA for initial positions and velocities.
- **Random new planets**: Added new planets with random characteristics.
- **Asteroids**: Added probability of appearance of randomly generated asteroids.
- **Collision detection**: The program detects collisions between bodies and merges them using the *Conservation of momentum*.

---

## Project Structure

- **Astrophysics_project_en.py**: Contains project code with english comments
- **Astrophysics_project_fr.py**: Contains project code with french comments
- **Project_Report_en**: Explains project parts,  architecture, laws studied in english
- **Rapport_Projet_fr**: Explains project parts,  architecture, laws studied in french
- **Stars.csv**: Dataset of stellar measurements.
- **Solar_System_data**: Physical and orbital properties of the Solar System planets.
- **Jupiter_Moons**: Orbital and physical characteristics of Jupiter's moons.
- **Exoplanets**: Catalog of exoplanets and their observed properties.

---

## Installation

1. Clone the repository:

```bash
git clone https://github.com/ali-hue-byte/Astrophysics_project
cd Astrophysics_project
```
2. Install the required dependencies:

```bash
python -m pip install matplotlib numpy csv random
```
3. Run the main application:

```bash
python Astrophysics_project_en.py
```

or 

```bash
python Astrophysics_project_fr.py
```
---
