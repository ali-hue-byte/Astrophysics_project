##### PROJET : Astrophysique et Astronomie : analyse et simulation #####


# Importation des modules
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import numpy as np
import csv

### Première partie : Analyse

# Fonction qui calcule la moyenne
def avg(x):
    return sum(x)/len(x)

# Fonction qui calcul l'écart-type
def ecart(x):
    variance = sum([(y-avg(x))**2 for y in x]) / len(x)
    return variance ** 0.5

# Fonction qui détecte les valeurs aberrantes
def aberr(x, y, a, b):
    residus = abs(y - (a * x + b))
    points_aberr = []

    moy_x = avg(x)
    moy_y = avg(y)

    sx = ecart(x)
    sy = ecart(y)

    for i in range(len(x)) :
        if (abs(x[i] - moy_x) > 2 * sx or abs(y[i]-moy_y) > 2 * sy) and residus[i] > 2 * ecart(residus):
            points_aberr.append(i) # renvoyer indice
    return points_aberr

# Fonction qui effectue la régression linéaire
def regression(x,y):
    mx = avg(x)
    my = avg(y)

    a = np.dot(x - mx, y - my) / np.dot(x - mx, x - mx)
    b = my - a * mx
    return a, b

## Loi de Kepler (lunes de jupiter) ##
# Source des données: https://gist.github.com/adrn/8191304#file-data-csv-L3

# Constantes
G = 6.67428e-11  # Constante de gravité en m³/(kg*s²)
G_u = G * ((1e-3)**3) * ((3600*24*365) ** 2) # Constante de gravité en km³/(kg*an²)
UA = 1.496e11 # 1UA en m
Mj = 1.898e27 # masse de jupiter en kg
Ms = 1.989e30 # masse de soleil en kg
an = 3.154e7 # 1an en s

# Lecture de fichier csv contenant les données
with open("Jupiter_Moons.csv","rt") as f:
    reader = csv.DictReader(f)
    data = [[x["Semi-major axis (km)"], x["Orbital period (days)"]] for x in reader]

x = np.array([float(i[0]) for i in data ]) # Demi-grand axe de l’orbite (km)
y = np.array([float(i[1])/365 for i in data ]) # Périod en années
d3 = x**3
T2 = y**2

a , b = regression(d3,T2)

print(f"## Loi de de Kepler (lunes je jupiter) ##\nCoefficient théorique: {(4*(np.pi**2))/(G_u*Mj)}")
print(f"Résultat de régression \na = {a}")
print(f"b = {b}")
# T² = a * d³
plt.scatter(d3,T2)
plt.plot(d3, a * d3 + b, label = f"Régression: T² = {a:.2e} * d³ + {b:.2e}", color = "black")
plt.title("Loi de Kepler (Lunes je jupiter)")
plt.xlabel("d³ [km³]")
plt.ylabel("T² [années²]")
plt.legend()
plt.show()


## Loi de Kepler (planètes du système solaire) ##
# Source des données: https://dataherb.github.io/flora/planets_in_solar_system/
with open("Solar_System_data.csv","rt") as f:
    reader = csv.DictReader(f)
    ls = [[r["semi_major_axis"], r["orbital_period"]] for r in reader]

d3_p = np.array([float(x[0])**3 for x in ls]) # d³ avec d en UA
T3_p = np.array([float(x[1])**2 for x in ls]) # T² avec T en années

a_p, b_p = regression(d3_p,T3_p)
# En utilisant les unités astronomiques (UA) et les années, la constante de la loi de kepler se simplifie en 1
print(f"\n## Loi de de Kepler (Planètes de système solaire) ##\nCoefficient théorique (en UA et années): {((4 * (np.pi**2)) / (G*Ms)) * (UA ** 3) / (an**2)}")
print(f"Résultat de régression \na = {a_p}")
print(f"b = {b_p}")

plt.scatter(d3_p,T3_p)
plt.plot(d3_p, a_p*d3_p+b_p, label = f"Régression: T² = {a_p:.2f} * d³ + {b_p:.2f}", color = "black")
plt.xlabel("d³ [UA³]")
plt.ylabel("T² [années²]")
plt.title("Loi de Kepler (Planètes de système solaire)")
plt.legend()
plt.show()

##

## Relation entre température et distance ##
# T est proportionnelle à d^-0.5
# Source des données: https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PS
with open("Exoplanets.csv", "rt") as f:
    data = csv.DictReader(f)
    ls_p = [[float(r["pl_orbsmax"]),  float(r["pl_eqt"]),float(r["pl_bmasse"])] for r in data if (r["pl_orbsmax"] != "" and r["pl_eqt"] != "" and r["pl_bmasse"] != "")]

# Séléction des exoplanètes massives (masse 5x plus grandes que le soleil)
# et à une distance < 4 UA
# pour éliminer les planètes gazeuses très éloignées de leur étoile et celles petites dont la température dépend fortement
# des facteurs locaux (albédo et luminosité)
x2 = np.array([i[0] for i in ls_p if i[0] < 4 and i[2] > 5]) # Distance entre les planètes et leurs étoiles (UA)
y2 = np.array([i[1] for i in ls_p if i[0] < 4 and i[2] > 5]) # Température (Kelvin)

# Transformation en échelle log10 pour linéariser la relation de puissance
log_x2 = np.log10(x2)
log_y2 = np.log10(y2)

# Régression linéaire: log(T) = a * log(D) + b
a2,b2 = regression(log_x2,log_y2)

# log(T) = a * log(D) + b  --> log(T) = log(D^a) + log(10^b) --> log(T) = log(10^b * D^a) --> T = 10^b * D^a
print("\n## Relation entre température et distance ##")
print(f"Relation obtenue par régression: T = {10**b2} * d^{a2}")

fig, axs = plt.subplots(2)
# Graphique 1: Echelle linéaire
axs[0].scatter(x2,y2) # Données réelles
x_dess = np.linspace(0.001,4, 1000)
axs[0].plot(x_dess,10**b2 * x_dess**a2, color="black", label = f"Régression: T = {10**b2:.2f} * d^({a2:.2f})") # Régression
# Graphique 2: Echelle logarithmique
axs[1].scatter(log_x2,log_y2) # Données réelles
axs[1].plot(log_x2, a2*log_x2 + b2, color="black", label = f"Régression: log(T) = {a2:.2f} * log(d) + {b2:.2f}") # Régression

axs[0].legend()
axs[1].legend()

axs[0].set_xlabel("Distance (d)")
axs[0].set_ylabel("Température (T)")

axs[1].set_xlabel("log10(d)")
axs[1].set_ylabel("log10(T)")

plt.suptitle("Relation entre température et distance")
plt.show()

##

## Relation entre masse et rayon (étoiles) ##
# Source des données: https://vizier.u-strasbg.fr/viz-bin/VizieR-4
with open("Stars.csv", "rt") as f:
    data = csv.DictReader(f)
    ls_p2 = [[float(r["Mass"]),  float(r["Rad"])] for r in data if (r["Mass"] != "" and r["Rad"] != "")]

# Etude des étoiles plus massives que le soleil
# Relation: R = M^0.57

# Séléction des étoiles dans la masse est supérieure à celle du soleil
x_ls_grand = np.array([i[0] for i in ls_p2 if i[0] > 1 ]) # Masse
y_ls_grand = np.array([i[1] for i in ls_p2 if i[0] > 1 ]) # Rayon

# Transformation en échelle log10 pour linéariser la relation de puissance
log_x_ls_grand = np.log10(x_ls_grand)
log_y_ls_grand = np.log10(y_ls_grand)

# Détection des valeurs aberrantes
indices_aberr = aberr(log_x_ls_grand,log_y_ls_grand, 0.57, 0)

# Enlever le points aberrants
n_log_x_ls_grand = np.array([log_x_ls_grand[i] for i in range(len(log_x_ls_grand)) if i not in indices_aberr])
n_log_y_ls_grand = np.array([log_y_ls_grand[i] for i in range(len(log_y_ls_grand)) if i not in indices_aberr])

n_x_ls_grand = np.array([x_ls_grand[i] for i in range(len(x_ls_grand)) if i not in indices_aberr])
n_y_ls_grand = np.array([y_ls_grand[i] for i in range(len(y_ls_grand)) if i not in indices_aberr])

# Régression linéaire: log(R) = a * log(M) + b
a3, b3 = regression(n_log_x_ls_grand,n_log_y_ls_grand)
print("\n## Relation entre masse et rayon des étoiles plus massives que le Soleil ##")
print(f"Relation obtenue par régression: R = {10**b3} * M^{a3}")

fig2, axs2 = plt.subplots(2)

# Graphique 1: Echelle linéaire
x_dess2 = np.linspace(1,30, 10000)
axs2[0].scatter(n_x_ls_grand,n_y_ls_grand)
axs2[0].plot(x_dess2,10**b3 * x_dess2**a3, color="black", label = f"Régression: R = {10**b3:.2f} * M^({a3:.2f})")

# Graphique 2: Echelle log10
axs2[1].scatter(n_log_x_ls_grand, n_log_y_ls_grand)
axs2[1].scatter([log_x_ls_grand[i] for i in indices_aberr], [log_y_ls_grand[i] for i in indices_aberr], color = "red", label = "Valeurs aberrantes")
axs2[1].plot(log_x_ls_grand, a3*log_x_ls_grand + b3, color="black", label = f"Régression: log(R) =  {a3:.2f} * log(M) + {b3:.2f}")

axs2[0].legend()
axs2[1].legend()

axs2[0].set_xlabel("Masse (M)")
axs2[0].set_ylabel("Rayon (R)")

axs2[1].set_xlabel("log10(M)")
axs2[1].set_ylabel("log10(R)")

plt.suptitle("Relation entre masse et rayon (étoiles plus massives que le soleil)")
plt.show()

# Etoiles moins massives que le soleil
# Relation: R = M^0.8
x_ls_petit = np.array([i[0] for i in ls_p2 if i[0] < 1 ]) # Masse
y_ls_petit = np.array([i[1] for i in ls_p2 if i[0] < 1 ]) # Rayon

# Transformation en échelle log10 pour linéariser la relation de puissance
log_x_ls_petit = np.log10(x_ls_petit)
log_y_ls_petit = np.log10(y_ls_petit)

# Détection des valeures aberrantes
indices = aberr(log_x_ls_petit,log_y_ls_petit, 0.8, 0)

# Enlever le points aberrants
log_x_ls_petit_filt = np.array([log_x_ls_petit[i] for i in range(len(log_x_ls_petit)) if i not in indices])
log_y_ls_petit_filt = np.array([log_y_ls_petit[i] for i in range(len(log_y_ls_petit)) if i not in indices])

x_ls_petit_filt = np.array([x_ls_petit[i] for i in range(len(x_ls_petit)) if i not in indices])
y_ls_petit_filt = np.array([y_ls_petit[i] for i in range(len(y_ls_petit)) if i not in indices])

# Régression linéaire: log(R) = a * log(M) + b
a4, b4 = regression(log_x_ls_petit_filt,log_y_ls_petit_filt)
print("\n## Relation entre masse et rayon des étoiles moins massives que le Soleil ##")
print(f"Relation obtenue par régression: R = {10**b4} * M^{a4}")

fig3, axs3 = plt.subplots(2)
# Graphique 1: Echelle linéaire
x_dess3 = np.linspace(0.001,1, 1000)
axs3[0].scatter(x_ls_petit_filt,y_ls_petit_filt)
axs3[0].plot(x_dess3,10**b4 * x_dess3**a4, color="black", label=f"Régression: R = {10**b4:.2f} * M^({a4:.2f})")

# Graphique 2: Echelle log10
axs3[1].scatter(log_x_ls_petit_filt, log_y_ls_petit_filt)
axs3[1].plot(log_x_ls_petit_filt, a4 * log_x_ls_petit_filt + b4, color="black", label = f"Régression: log(R) =  {a4:.2f} * log(M) + {b4:.2f}")
axs3[1].scatter([log_x_ls_petit[i] for i in indices], [log_y_ls_petit[i] for i in indices], color = "red", label = "Valeurs aberrantes")

axs3[0].legend()
axs3[1].legend()

axs3[0].set_xlabel("Masse (M)")
axs3[0].set_ylabel("Rayon (R)")

axs3[1].set_xlabel("log10(M)")
axs3[1].set_ylabel("log10(R)")

plt.suptitle("Relation entre masse et rayon (étoiles moins massives que le soleil)")
plt.show()
##

### Deuxième partie : simulation
print("\n##### SIMULATION #####\n")

dt = 3600*24 # 1 jour

# date de début du simulation
hour = 0
day = 1
month = 1
year = 2026

# Données réelles de 1 janvier 2026 (https://ssd.jpl.nasa.gov/horizons/app.html#/)
planets = [
    {
        "nom": "Soleil",
        "r": 695700 * 1000,
        "m": 1.988410e30,
        "p": np.array([[0, 0]]),
        "v": np.array([0, 0]),
        "color":"yellow"
    },
    {
        "nom": "Mercure",
        "r": 2439.7e3,
        "m": 3.302e23,
        "p": np.array([[-3.265252087416521e+07 * 1000, -6.204436222997473e+07 * 1000]]),
        "v": np.array([3.331154717745245e+01 * 1000, -2.032288538992848e+01 * 1000]),
        "color":"dimgray"
    },
    {
        "nom": "Venus",
        "r": 6051.8e3,
        "m": 48.685e23,
        "p": np.array([[1.283698332887579e+07 * 1000, -1.088018897929847e+08 * 1000]]),
        "v": np.array([3.453518141145980e+01 * 1000, 4.156437812065552E+00 * 1000]),
        "color":"gold"
    },
    {
        "nom": "Terre",
        "r": 6371.0e3,
        "m": 5.97219e24,
        "p": np.array([[-2.653100241556548e+07 * 1000, 1.439468995740296e+08 * 1000]]),
        "v": np.array([-2.977650610770464e+01 * 1000, -5.395962660572101e+00 * 1000]),
        "color":"dodgerblue"
    },
    {
        "nom": "Mars",
        "r": 3389.5e3,
        "m": 6.4171e23,
        "p": np.array([[5.049113049487789e+07 * 1000, -2.083203224890075e+08 * 1000]]),
        "v": np.array([2.445957741742463e+01 * 1000, 7.861473739133970e+00 * 1000]),
        "color":"red"
    },
    {
        "nom": "Jupiter",
        "r": 69911e3,
        "m": 18.9819e26,
        "p": np.array([[-2.538782093363155e+08 * 1000, 7.365225315477104e+08 * 1000]]),
        "v": np.array([-1.250761633714600E+01 * 1000, -3.639986644887777e+00 * 1000]),
        "color":"peru"
    },
    {
        "nom": "Saturne",
        "r": 58232e3,
        "m": 5.6834e26,
        "p": np.array([[1.421819203232436e+09 * 1000, 3.772943733460353e+07 * 1000]]),
        "v": np.array([-7.898651748040556e-01 * 1000, 9.634221220286250e+00 * 1000]),
        "color":"wheat"
    },
    {
        "nom": "Uranus",
        "r": 25362e3,
        "m": 86.813e24,
        "p": np.array([[1.477614612946379e+09 * 1000, 2.512418267148477e+09 * 1000]]),
        "v": np.array([-5.920365479337566e+00 * 1000, 3.134957028375065e+00 * 1000]),
        "color":"cyan"
    },
    {
        "nom": "Neptune",
        "r": 24622e3,
        "m": 102.409e24,
        "p": np.array([[4.468346746922093e+09 * 1000, 7.680447071835697e+07 * 1000]]),
        "v": np.array([-1.294614448856751e-01 * 1000, 5.465955556668026e+00 * 1000]),
        "color":"navy"
    },
]

# Fonction pour générer des objets aléatoirement
def generer_aleatoire(choice):
    global planets

    if choice == "Planète":
        minimum = 1e23   # masse minimale
        maximum = 1e27   # masse maximale
        color = "green"
        densite = random.randrange(1000, 6000) # kg/m^3
    else:
        minimum = 1e12
        maximum = 1e18
        color = "gray"
        densite = random.randrange(2000, 9000)  # kg/m^3

    nom = "Aleatoire " + str(len(planets) - 9)
    masse = random.randrange(int(minimum),int(maximum)) # masse aléatoire selon type d'objet
    distan = random.randrange(int(1e10),int(5e12)) # distance du soleil aléatoire
    angle = 2 * np.pi * random.random() # Angle aléatoire (entre 0 et 2pi) pour une répartition réaliste
    angle += random.uniform(-2, 2) # (float aléatoire) petite perturbation dans l'angle
    position = [distan * np.cos(angle), distan * np.sin(angle)] # Calcul des cordonnées x et y de façon
                                                                # à suivre un orbit un peu réaliste
    vitesse = np.sqrt(G * planets[0]["m"] / distan) * random.uniform(0.8,1.2) # Norme de la vitesse à partir de distance et masse
                                                    # v = racine(GM/d)
    direction = random.choice([-1, 1])
    vx = direction * (-distan * np.sin(angle) / distan) * vitesse  # vx = dx/dt et vy = dy/dt (x = dist*cos(angle) et y = dist*sin(angle) )
    vy = direction * (distan * np.cos(angle) / distan) * vitesse   # on divise par la distance pour obtenir un vecteur de norme 1
                                            # puis, on multiplie par la norme de la vitesse et direction aléatoire

    if choice == "Astéroide":                        # Ajouter des déformations pour les astéroïdes
        position[0] += random.uniform(-1e8, 1e8)
        position[1] += random.uniform(-1e8, 1e8)
        vx *= random.uniform(0.8,1.2)
        vy *= random.uniform(0.8,1.2)


    # Ajout du nouvel objet
    planets.append({
        "nom": nom,
        "r": (3 * masse / (4 * np.pi * densite)) ** (1 / 3),   # m = densité * 4/3 * pi * r^3
        "m": masse,
        "p": np.array([position]),
        "v": np.array([vx, vy]),
        "color":color
    })

    # Affichage des caractéristiques du nouvel objet
    print(f"{choice} généré(e): {nom}\nMasse: {masse:.2e}\nPosition: {position}\nVitesse: {[vx,vy]}")
    print("######################################")

# Générer des planètes aléatoires
for i in range(5):
   generer_aleatoire("Planète")

# Fonction qui calcule la distance entre objets i et j
def distance(i,j):
    x = planets[i]["p"][-1][0] - planets[j]["p"][-1][0]
    y = planets[i]["p"][-1][1] - planets[j]["p"][-1][1]

    return np.sqrt(x**2 + y**2)

# Fonction qui calcule accélération d'un objet
# aj = 1/mj * somme des forces
# vect(force) = G*(mi*mj)/r^2 * vect(u)
# avec vect(u) = vecteur unitaire = vect(r)/r
# vect(force) = G*(mi*mj)*vect(r) / r^3
# et donc a = G * m_i * vect(r) / r^3

def acceleration(j):
    a = np.array([0.0, 0.0])
    for i in range(len(planets)):
        if i != j:
            ai = G * planets[i]["m"] / distance(i, j) ** 3 * (planets[i]["p"][-1] - planets[j]["p"][-1])
            a += ai
    return a

# Création de la figure

plt.style.use('dark_background')
fig = plt.figure(figsize=(15, 8))
lines = [plt.plot([], [])[0] for i in range(len(planets))] # stocke les trajectoires (historique des positions)
points = [plt.scatter([], []) for i in range(len(planets))] # points : la position actuelle de chaque corps
plt.xlim(-5e12, 5e12)
plt.ylim(-5e12, 5e12)
date = plt.text(0.02, 0.95, "Date: ", transform=plt.gca().transAxes)
plt.title("Simulation du mouvement des planètes")

# Fonction pour vérifier qu'une année est bissextile
def leap(i):
    if i % 400 == 0:
        return True
    elif i % 100 == 0:
        return False
    elif i % 4 == 0:
        return True
    else:
        return False

# Fonction qui fusionne deux objets
def fusion(surv, disparait):

    m1 = float(planets[surv]["m"])
    m2 = float(planets[disparait]["m"])

    v1 = np.array(planets[surv]["v"])
    v2 = np.array(planets[disparait]["v"])

    r1 = planets[surv]["r"]
    r2 = planets[disparait]["r"]

    m_new = m1 + m2
    v_new = (m1 * v1 + m2 * v2) / m_new # conservation de quantité de mouvement m1*v1 + m2*v2 = (m1+m2)v

    planets[surv]["m"] = m_new
    planets[surv]["v"] = v_new

    # approximation : conservation du volume total
    # V = 4/3 * pi * r^3
    # Vt = V1 + V2 nous donne que rt^3 = r1^3 + r2^3
    planets[surv]["r"] = (r1 ** 3 + r2 ** 3) ** (1 / 3)


def animate(i):
    p = 0.01   # probabilité d'apparition d'un astéroide

    global planets, day, month, year, lines, points, hour

    if random.random() < p:
        generer_aleatoire("Astéroide") # générer astéroide aléatoire
        lines.append(plt.plot([], [])[0]) # Ajout du graphe de l'astéroide
        points.append(plt.scatter([], []))

    # Méthode de Verlet

    a = np.array([acceleration(j) for j in range(len(planets))]) # calcul d'accélération à l'instant t
    for n in range(len(planets)):
        pos = a[n]/2 * dt**2 + planets[n]["v"]*dt + planets[n]["p"][-1] # Calcul de position à l'instant t + dt
        planets[n]["p"] = np.append(planets[n]["p"], [pos], axis=0)  # x(t+dt) = 1/2*a(t)*t² + v(t)*t + x(t)

    a_new = np.array([acceleration(j) for j in range(len(planets))]) # calcul de l'accélération à l'instant t + dt


    for k in range(len(planets)):
        vitesse = planets[k]["v"] + ((a[k] + a_new[k]) * dt) / 2 # Calcul de vitesses à l'instant t + dt (dépend de la vitesse à
                                                                 # l'instant t et l'accélération à l'instant t et t+dt)
        planets[k]["v"] = vitesse

    # Vérifier si il y'a des collisions
    planetes_a_disparaitre = []
    for k in range(len(planets)):
        for j in range(k+1, len(planets)):
            dist = planets[k]["r"] + planets[j]["r"]

            if distance(k,j) < dist :
                print("\n######################################################")
                print(f"Collision entre {planets[k]['nom']} et {planets[j]['nom']}")
                print("######################################################\n")
                if planets[k]["m"] > planets[j]["m"]:
                    surv = k
                    disparait = j
                else:
                    surv = j
                    disparait = k
                fusion(surv, disparait) # fusionner les deux objets
                planetes_a_disparaitre.append(disparait)

    asteroide = False
    Planete = False
    # Mise à jour des graphes
    for j in range(len(planets)):

        u = planets[j]["p"][:, 0]
        v = planets[j]["p"][:, 1]
        if j < len(lines) :  # Eviter les erreurs
            lines[j].set_data(u, v)
            points[j].set_offsets([planets[j]["p"][-1][0], planets[j]["p"][-1][1]])

            # Taille des points dépend du rayon (sauf le soleil car il va apparaitre trop gros)
            points[j].set_sizes([planets[j]["r"] / 2e5]) if j != 0 else points[j].set_sizes([40])
            lines[j].set_color(planets[j]["color"])
            points[j].set_color(planets[j]["color"])


        if j == 0 :
            points[j].set_label(planets[j]["nom"])
        elif j > 8 :
            if (planets[j]["m"] >= 1e20) and (not Planete):
                Planete = True
                lines[j].set_label("Planète")
            elif (not asteroide) and (planets[j]["m"] <= 1e18):
                asteroide = True
                lines[j].set_label("Astéroide")
        else :
            lines[j].set_label(planets[j]["nom"])

        plt.legend()


    # Date actuelle
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




    # Enlever les positions très anciennes après 500 pas
    for k in planets:
         if len(k["p"]) > 500:
            k["p"] = k["p"][1:]
    date.set_text(f"Date: {int(day)}/{int(month)}/{int(year)}") # changer la date

    # Enlever les objets collisionées et laisser que la fusion
    for j in sorted(set(planetes_a_disparaitre), reverse=True):
        lines[j].remove()
        points[j].remove()
        lines.pop(j)
        points.pop(j)
        planets.pop(j)

    return lines

ani = animation.FuncAnimation(fig, animate,interval=10, blit=False, repeat=False)
plt.show()