##### PROJECT : PROJECT: Astrophysics and Astronomy: Analysis and Simulation #####


# Module imports
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import numpy as np
import csv

### Part 1: Analysis

# Computes the mean of a dataset
def avg(x):
    return sum(x)/len(x)

# Computes the standard deviation of a dataset
def ecart(x):
    variance = sum([(y-avg(x))**2 for y in x]) / len(x)
    return variance ** 0.5

# Detects outliers based on residuals and standard deviation
def aberr(x, y, a, b):
    residus = abs(y - (a * x + b))
    points_aberr = []

    moy_x = avg(x)
    moy_y = avg(y)

    sx = ecart(x)
    sy = ecart(y)

    for i in range(len(x)) :
        if (abs(x[i] - moy_x) > 2 * sx or abs(y[i]-moy_y) > 2 * sy) and residus[i] > 2 * ecart(residus):
            points_aberr.append(i) # returns index
    return points_aberr

# Performs linear regression using least squares method
def regression(x,y):
    mx = avg(x)
    my = avg(y)

    a = np.dot(x - mx, y - my) / np.dot(x - mx, x - mx)
    b = my - a * mx
    return a, b

## Kepler's Law (Jupiter's moons) ##
# Data source: https://gist.github.com/adrn/8191304#file-data-csv-L3

# Constants
G = 6.67428e-11  # Gravitational constant in m³/(kg*s²)
G_u = G * ((1e-3)**3) * ((3600*24*365) ** 2) # Gravitational constant in km³/(kg*year²)
UA = 1.496e11 # 1 AU in m
Mj = 1.898e27 # Jupiter mass in kg
Ms = 1.989e30 # Sun mass in kg
an = 3.154e7  # 1 year in s

# Read CSV file containing orbital data
with open("Jupiter_Moons.csv","rt") as f:
    reader = csv.DictReader(f)
    data = [[x["Semi-major axis (km)"], x["Orbital period (days)"]] for x in reader]

x = np.array([float(i[0]) for i in data ]) # Semi-major axis of the orbit (km)
y = np.array([float(i[1])/365 for i in data ]) # Orbital period in years
d3 = x**3
T2 = y**2

a , b = regression(d3,T2)

print(f"## Kepler's Law (Jupiter's moons) ##\nTheoretical coefficient: {(4*(np.pi**2))/(G_u*Mj)}")
print(f"Regression result \na = {a}")
print(f"b = {b}")
# T² = a * d³
plt.scatter(d3,T2)
plt.plot(d3, a * d3 + b, label = f"Regression: T² = {a:.2e} * d³ + {b:.2e}", color = "black")
plt.title("Kepler's Law (Jupiter's moons)")
plt.xlabel("d³ [km³]")
plt.ylabel("T² [years²]")
plt.legend()
plt.show()


## Kepler's law (Solar System planets) ##
# Data source: https://dataherb.github.io/flora/planets_in_solar_system/
with open("Solar_System_data.csv","rt") as f:
    reader = csv.DictReader(f)
    ls = [[r["semi_major_axis"], r["orbital_period"]] for r in reader]

d3_p = np.array([float(x[0])**3 for x in ls]) # d³ with d in AU
T3_p = np.array([float(x[1])**2 for x in ls]) # T² with T in years

a_p, b_p = regression(d3_p,T3_p)
# Using astronomical units (AU) and years, Kepler's constant simplifies to 1
print(f"\n## Kepler's Law (Solar System planets) ##\nTheoretical coefficient (in AU and years): {((4 * (np.pi**2)) / (G*Ms)) * (UA ** 3) / (an**2)}")
print(f"Regression result \na = {a_p}")
print(f"b = {b_p}")

plt.scatter(d3_p,T3_p)
plt.plot(d3_p, a_p*d3_p+b_p, label = f"Regression: T² = {a_p:.2f} * d³ + {b_p:.2f}", color = "black")
plt.xlabel("d³ [AU³]")
plt.ylabel("T² [years²]")
plt.title("Kepler's Law (Solar System planets)")
plt.legend()
plt.show()

##

## Temperature-Distance Relationship ##
# T is proportional to d^-0.5
# Data source: https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PS
with open("Exoplanets.csv", "rt") as f:
    data = csv.DictReader(f)
    ls_p = [[float(r["pl_orbsmax"]),  float(r["pl_eqt"]),float(r["pl_bmasse"])] for r in data if (r["pl_orbsmax"] != "" and r["pl_eqt"] != "" and r["pl_bmasse"] != "")]

# Filter: keep only massive exoplanets (mass > 5 Earth masses) within 4 AU
# to filter out gas giants far from their star and small planets whose
# temperature strongly depends on local factors (albedo and luminosity)
x2 = np.array([i[0] for i in ls_p if i[0] < 4 and i[2] > 5]) # Distance (AU)
y2 = np.array([i[1] for i in ls_p if i[0] < 4 and i[2] > 5]) # Temperature (Kelvin)

# Log10 transformation to linearize the power-law relationship
log_x2 = np.log10(x2)
log_y2 = np.log10(y2)

# Linear regression: log(T) = a * log(D) + b
a2,b2 = regression(log_x2,log_y2)

# log(T) = a * log(D) + b  --> log(T) = log(D^a) + log(10^b) --> log(T) = log(10^b * D^a) --> T = 10^b * D^a
print("\n## Temperature-Distance Relationship ##")
print(f"Regression result: T = {10**b2} * d^{a2}")

fig, axs = plt.subplots(2)
# Graph 1: Linear scale
axs[0].scatter(x2,y2) # Real data
x_dess = np.linspace(0.001,4, 1000)
axs[0].plot(x_dess,10**b2 * x_dess**a2, color="black", label = f"Regression: T = {10**b2:.2f} * d^({a2:.2f})") # Regression
# Graph 2: Logarithmic scale
axs[1].scatter(log_x2,log_y2) # Real data
axs[1].plot(log_x2, a2*log_x2 + b2, color="black", label = f"Regression: log(T) = {a2:.2f} * log(d) + {b2:.2f}") # Regression

axs[0].legend()
axs[1].legend()

axs[0].set_xlabel("Distance (d)")
axs[0].set_ylabel("Temperature (T)")

axs[1].set_xlabel("log10(d)")
axs[1].set_ylabel("log10(T)")

plt.suptitle("Temperature-Distance Relationship")
plt.show()

##

## Mass-Radius Relationship (stars) ##
# Data source: https://vizier.u-strasbg.fr/viz-bin/VizieR-4
with open("Stars.csv", "rt") as f:
    data = csv.DictReader(f)
    ls_p2 = [[float(r["Mass"]),  float(r["Rad"])] for r in data if (r["Mass"] != "" and r["Rad"] != "")]

# Study of stars more massive than the Sun
# R = M^0.57

# Select stars whose mass is greater than the Sun's
x_ls_grand = np.array([i[0] for i in ls_p2 if i[0] > 1 ]) # Mass
y_ls_grand = np.array([i[1] for i in ls_p2 if i[0] > 1 ]) # Radius

# Transformation to logarithmic scale to linearize the power-law relationship
log_x_ls_grand = np.log10(x_ls_grand)
log_y_ls_grand = np.log10(y_ls_grand)

# Outlier detection
indices_aberr = aberr(log_x_ls_grand,log_y_ls_grand, 0.57, 0)

# Remove outlier points
n_log_x_ls_grand = np.array([log_x_ls_grand[i] for i in range(len(log_x_ls_grand)) if i not in indices_aberr])
n_log_y_ls_grand = np.array([log_y_ls_grand[i] for i in range(len(log_y_ls_grand)) if i not in indices_aberr])

n_x_ls_grand = np.array([x_ls_grand[i] for i in range(len(x_ls_grand)) if i not in indices_aberr])
n_y_ls_grand = np.array([y_ls_grand[i] for i in range(len(y_ls_grand)) if i not in indices_aberr])

# Linear regression: log(R) = a * log(M) + b
a3, b3 = regression(n_log_x_ls_grand,n_log_y_ls_grand)
print("\n## Mass-Radius Relationship (stars more massive than the Sun) ##")
print(f"Regression result: R = {10**b3} * M^{a3}")

fig2, axs2 = plt.subplots(2)

# Graph 1: Linear scale
x_dess2 = np.linspace(1,30, 10000)
axs2[0].scatter(n_x_ls_grand,n_y_ls_grand)
axs2[0].plot(x_dess2,10**b3 * x_dess2**a3, color="black", label = f"Regression: R = {10**b3:.2f} * M^({a3:.2f})")

# Graph 2: Logarithmic scale
axs2[1].scatter(n_log_x_ls_grand, n_log_y_ls_grand)
axs2[1].scatter([log_x_ls_grand[i] for i in indices_aberr], [log_y_ls_grand[i] for i in indices_aberr], color = "red", label = "Outliers")
axs2[1].plot(log_x_ls_grand, a3*log_x_ls_grand + b3, color="black", label = f"Regression: log(R) = {a3:.2f} * log(M) + {b3:.2f}")

axs2[0].legend()
axs2[1].legend()

axs2[0].set_xlabel("Mass (M)")
axs2[0].set_ylabel("Radius (R)")

axs2[1].set_xlabel("log10(M)")
axs2[1].set_ylabel("log10(R)")

plt.suptitle("Mass-Radius Relationship (stars more massive than the Sun)")
plt.show()

# Stars less massive than the Sun
# R = M^0.8
x_ls_petit = np.array([i[0] for i in ls_p2 if i[0] < 1 ]) # Mass
y_ls_petit = np.array([i[1] for i in ls_p2 if i[0] < 1 ]) # Radius

# Transformation to logarithmic scale to linearize the power-law relationship
log_x_ls_petit = np.log10(x_ls_petit)
log_y_ls_petit = np.log10(y_ls_petit)

# Outlier detection
indices = aberr(log_x_ls_petit,log_y_ls_petit, 0.8, 0)

# Remove outlier points
log_x_ls_petit_filt = np.array([log_x_ls_petit[i] for i in range(len(log_x_ls_petit)) if i not in indices])
log_y_ls_petit_filt = np.array([log_y_ls_petit[i] for i in range(len(log_y_ls_petit)) if i not in indices])

x_ls_petit_filt = np.array([x_ls_petit[i] for i in range(len(x_ls_petit)) if i not in indices])
y_ls_petit_filt = np.array([y_ls_petit[i] for i in range(len(y_ls_petit)) if i not in indices])

# Linear regression: log(R) = a * log(M) + b
a4, b4 = regression(log_x_ls_petit_filt,log_y_ls_petit_filt)
print("\n## Mass-Radius Relationship (stars less massive than the Sun) ##")
print(f"Regression result: R = {10**b4} * M^{a4}")

fig3, axs3 = plt.subplots(2)
# Graph 1: Linear scale
x_dess3 = np.linspace(0.001,1, 1000)
axs3[0].scatter(x_ls_petit_filt,y_ls_petit_filt)
axs3[0].plot(x_dess3,10**b4 * x_dess3**a4, color="black", label=f"Regression: R = {10**b4:.2f} * M^({a4:.2f})")

# Graph 2: Logarithmic scale
axs3[1].scatter(log_x_ls_petit_filt, log_y_ls_petit_filt)
axs3[1].plot(log_x_ls_petit_filt, a4 * log_x_ls_petit_filt + b4, color="black", label = f"Regression: log(R) = {a4:.2f} * log(M) + {b4:.2f}")
axs3[1].scatter([log_x_ls_petit[i] for i in indices], [log_y_ls_petit[i] for i in indices], color = "red", label = "Outliers")

axs3[0].legend()
axs3[1].legend()

axs3[0].set_xlabel("Mass (M)")
axs3[0].set_ylabel("Radius (R)")

axs3[1].set_xlabel("log10(M)")
axs3[1].set_ylabel("log10(R)")

plt.suptitle("Mass-Radius Relationship (stars less massive than the Sun)")
plt.show()
##

### Part 2 : Simulation
print("\n##### SIMULATION #####\n")

dt = 3600*24 # 1 day

# Simulation start date
hour = 0
day = 1
month = 1
year = 2026

# Real data from January 1st, 2026 (https://ssd.jpl.nasa.gov/horizons/app.html#/)
planets = [
    {
        "name": "Sun",
        "r": 695700 * 1000,
        "m": 1.988410e30,
        "p": np.array([[0, 0]]),
        "v": np.array([0, 0]),
        "color":"yellow"
    },
    {
        "name": "Mercury",
        "r": 2439.7e3,
        "m": 3.302e23,
        "p": np.array([[-3.265252087416521e+07 * 1000, -6.204436222997473e+07 * 1000]]),
        "v": np.array([3.331154717745245e+01 * 1000, -2.032288538992848e+01 * 1000]),
        "color":"dimgray"
    },
    {
        "name": "Venus",
        "r": 6051.8e3,
        "m": 48.685e23,
        "p": np.array([[1.283698332887579e+07 * 1000, -1.088018897929847e+08 * 1000]]),
        "v": np.array([3.453518141145980e+01 * 1000, 4.156437812065552E+00 * 1000]),
        "color":"gold"
    },
    {
        "name": "Earth",
        "r": 6371.0e3,
        "m": 5.97219e24,
        "p": np.array([[-2.653100241556548e+07 * 1000, 1.439468995740296e+08 * 1000]]),
        "v": np.array([-2.977650610770464e+01 * 1000, -5.395962660572101e+00 * 1000]),
        "color": "dodgerblue"
    },
    {
        "name": "Mars",
        "r": 3389.5e3,
        "m": 6.4171e23,
        "p": np.array([[5.049113049487789e+07 * 1000, -2.083203224890075e+08 * 1000]]),
        "v": np.array([2.445957741742463e+01 * 1000, 7.861473739133970e+00 * 1000]),
        "color": "red"
    },
    {
        "name": "Jupiter",
        "r": 69911e3,
        "m": 18.9819e26,
        "p": np.array([[-2.538782093363155e+08 * 1000, 7.365225315477104e+08 * 1000]]),
        "v": np.array([-1.250761633714600E+01 * 1000, -3.639986644887777e+00 * 1000]),
        "color": "peru"
    },
    {
        "name": "Saturn",
        "r": 58232e3,
        "m": 5.6834e26,
        "p": np.array([[1.421819203232436e+09 * 1000, 3.772943733460353e+07 * 1000]]),
        "v": np.array([-7.898651748040556e-01 * 1000, 9.634221220286250e+00 * 1000]),
        "color": "wheat"
    },
    {
        "name": "Uranus",
        "r": 25362e3,
        "m": 86.813e24,
        "p": np.array([[1.477614612946379e+09 * 1000, 2.512418267148477e+09 * 1000]]),
        "v": np.array([-5.920365479337566e+00 * 1000, 3.134957028375065e+00 * 1000]),
        "color": "cyan"
    },
    {
        "name": "Neptune",
        "r": 24622e3,
        "m": 102.409e24,
        "p": np.array([[4.468346746922093e+09 * 1000, 7.680447071835697e+07 * 1000]]),
        "v": np.array([-1.294614448856751e-01 * 1000, 5.465955556668026e+00 * 1000]),
        "color": "navy"
    },
]


# Generates a random planet or asteroid with realistic parameters
def generate_random(choice):
    global planets

    if choice == "Planet":
        minimum = 1e23  # minimum mass
        maximum = 1e27  # maximum mass
        color = "green"
        density = random.randrange(1000, 6000)  # kg/m^3
    else:
        minimum = 1e12
        maximum = 1e18
        color = "gray"
        density = random.randrange(2000, 9000)  # kg/m^3

    name = "Random " + str(len(planets) - 9)
    mass = random.randrange(int(minimum), int(maximum))  # random mass based on object type
    dist = random.randrange(int(1e10), int(5e12))  # random distance from the Sun
    angle = 2 * np.pi * random.random()  # random angle for realistic distribution
    angle += random.uniform(-2, 2)  # small perturbation in angle
    position = [dist * np.cos(angle), dist * np.sin(angle)]  # compute x and y coordinates for a realistic orbit
    speed = np.sqrt(G * planets[0]["m"] / dist) * random.uniform(0.8,
                                                                 1.2)  # orbital speed magnitude from distance and mass: v = sqrt(GM/d)
    direction = random.choice([-1, 1])
    vx = direction * (-dist * np.sin(
        angle) / dist) * speed  # vx = dx/dt and vy = dy/dt (x = dist*cos(angle) and y = dist*sin(angle))
    vy = direction * (dist * np.cos(angle) / dist) * speed  # divide by distance to get a unit vector (norm = 1)
    # then multiply by the speed magnitude and random direction

    if choice == "Asteroid":  # Add slight deformations for asteroids
        position[0] += random.uniform(-1e8, 1e8)
        position[1] += random.uniform(-1e8, 1e8)
        vx *= random.uniform(0.8, 1.2)
        vy *= random.uniform(0.8, 1.2)

    # Add the new object to the planets list
    planets.append({
        "name": name,
        "r": (3 * mass / (4 * np.pi * density)) ** (1 / 3),  # m = density * 4/3 * pi * r^3
        "m": mass,
        "p": np.array([position]),
        "v": np.array([vx, vy]),
        "color": color
    })

    # Print characteristics of the new object
    print(f"{choice} generated: {name}\nMass: {mass:.2e}\nPosition: {position}\nVelocity: {[vx, vy]}")
    print("######################################")


# Generate random planets
for i in range(5):
    generate_random("Planet")


# Computes the Euclidean distance between objects i and j
def distance(i, j):
    x = planets[i]["p"][-1][0] - planets[j]["p"][-1][0]
    y = planets[i]["p"][-1][1] - planets[j]["p"][-1][1]

    return np.sqrt(x ** 2 + y ** 2)


# Computes the gravitational acceleration of object j
# aj = 1/mj * sum of forces
# force_vector = G*(mi*mj)/r^2 * unit_vector
# with unit_vector = r_vector/r
# force_vector = G*(mi*mj)*r_vector / r^3
# therefore a = G * m_i * r_vector / r^3

def acceleration(j):
    a = np.array([0.0, 0.0])
    for i in range(len(planets)):
        if i != j:
            ai = G * planets[i]["m"] / distance(i, j) ** 3 * (planets[i]["p"][-1] - planets[j]["p"][-1])
            a += ai
    return a


# Figure setup

plt.style.use('dark_background')
fig = plt.figure(figsize=(15, 8))
lines = [plt.plot([], [])[0] for i in range(len(planets))]  # stores trajectories (position history)
points = [plt.scatter([], []) for i in range(len(planets))]  # points: current position of each body
plt.xlim(-5e12, 5e12)
plt.ylim(-5e12, 5e12)
date = plt.text(0.02, 0.95, "Date: ", transform=plt.gca().transAxes)
plt.title("Simulation of planetary motion")


# Checks if a year is a leap year
def leap(i):
    if i % 400 == 0:
        return True
    elif i % 100 == 0:
        return False
    elif i % 4 == 0:
        return True
    else:
        return False


# Merges two colliding objects using momentum conservation
def merge(survivor, vanishing):
    m1 = float(planets[survivor]["m"])
    m2 = float(planets[vanishing]["m"])

    v1 = np.array(planets[survivor]["v"])
    v2 = np.array(planets[vanishing]["v"])

    r1 = planets[survivor]["r"]
    r2 = planets[vanishing]["r"]

    m_new = m1 + m2
    v_new = (m1 * v1 + m2 * v2) / m_new  # conservation of momentum: m1*v1 + m2*v2 = (m1+m2)*v_new

    planets[survivor]["m"] = m_new
    planets[survivor]["v"] = v_new

    # → approximation: total volume conservation
    # → Vt = V1 + V2 → rt³ = r1³ + r2³

    planets[survivor]["r"] = (r1 ** 3 + r2 ** 3) ** (1 / 3)


def animate(i):
    p = 0.01  # probability of an asteroid appearing

    global planets, day, month, year, lines, points, hour

    if random.random() < p:
        generate_random("Asteroid")  # generate a random asteroid
        lines.append(plt.plot([], [])[0])
        points.append(plt.scatter([], []))

    # Verlet integration method

    a = np.array([acceleration(j) for j in range(len(planets))])  # compute acceleration at time t
    for n in range(len(planets)):
        pos = a[n] / 2 * dt ** 2 + planets[n]["v"] * dt + planets[n]["p"][-1]  # compute position at time t + dt
        planets[n]["p"] = np.append(planets[n]["p"], [pos], axis=0)  # x(t+dt) = 1/2*a(t)*t² + v(t)*t + x(t)

    a_new = np.array([acceleration(j) for j in range(len(planets))])  # compute acceleration at time t + dt

    for k in range(len(planets)):
        velocity = planets[k]["v"] + ((a[k] + a_new[k]) * dt) / 2  # compute velocity at time t + dt
        planets[k]["v"] = velocity

    # Check for collisions between all pairs
    planets_to_remove = []
    for k in range(len(planets)):
        for j in range(k + 1, len(planets)):
            dist = planets[k]["r"] + planets[j]["r"]

            if distance(k, j) < dist:
                print("\n######################################################")
                print(f"Collision between {planets[k]['name']} and {planets[j]['name']}")
                print("######################################################\n")
                if planets[k]["m"] > planets[j]["m"]:
                    survivor = k
                    vanishing = j
                else:
                    survivor = j
                    vanishing = k
                merge(survivor, vanishing)  # merge the two objects
                planets_to_remove.append(vanishing)

    asteroid = False
    planet = False
    # Update plots
    for j in range(len(planets)):

        u = planets[j]["p"][:, 0]
        v = planets[j]["p"][:, 1]
        if j < len(lines):
            lines[j].set_data(u, v)
            points[j].set_offsets([planets[j]["p"][-1][0], planets[j]["p"][-1][1]])

            # Point size depends on radius (except Sun)
            points[j].set_sizes([planets[j]["r"] / 2e5]) if j != 0 else points[j].set_sizes([40])
            lines[j].set_color(planets[j]["color"])
            points[j].set_color(planets[j]["color"])

        if j == 0:
            points[j].set_label(planets[j]["name"])
        elif j > 8:
            if (planets[j]["m"] >= 1e20) and (not planet):
                planet = True
                lines[j].set_label("Planet")
            elif (not asteroid) and (planets[j]["m"] <= 1e18):
                asteroid = True
                lines[j].set_label("Asteroid")
        else:
            lines[j].set_label(planets[j]["name"])

        plt.legend()

    # Current simulation date
    day += dt // (3600 * 24)
    hour += (dt % (3600 * 24)) / 3600

    if hour >= 24:
        day += int(hour // 24)
        hour = hour % 24

    if month == 2 and not leap(year) and day > 28:
        day -= 28
        month += 1
    elif month == 2 and leap(year) and day > 29:
        day -= 29
        month += 1
    elif month == 12 and day > 31:
        month = 1
        day -= 31
        year += 1
    elif day > 30 and month in [4, 6, 9, 11]:
        day -= 30
        month += 1
    elif day > 31 and month in [1, 3, 5, 7, 8, 10]:
        day -= 31
        month += 1

    # Remove old positions after 500 steps to save memory
    for k in planets:
        if len(k["p"]) > 500:
            k["p"] = k["p"][1:]
    date.set_text(f"Date: {int(day)}/{int(month)}/{int(year)}")  # update the date

    # Remove collided objects and keep only the merged result
    for j in sorted(set(planets_to_remove), reverse=True):
        lines[j].remove()
        points[j].remove()
        lines.pop(j)
        points.pop(j)
        planets.pop(j)

    return lines


ani = animation.FuncAnimation(fig, animate, interval=10, blit=False, repeat=False)
plt.show()