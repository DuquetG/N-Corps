import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Lire les données du fichier .dat
with open('bodies_movement.dat', 'r') as file:
    lines = file.readlines()

# Séparer les valeurs x et y pour chaque objet
data = np.array([list(map(float, line.split())) for line in lines])
print(data)

# Nombre d'objets
num_objects = len(data[0]) // 2

# Créer une figure et un axe
fig, ax = plt.subplots()
ax.set_xlim(-200, 200)
ax.set_ylim(-20000, 20000)

# Créer des objets pour chaque point initial
points, = ax.plot([], [], marker='o', linestyle='', color='b', markersize=10)
trails, = ax.plot([], [], linestyle='-', color='gray', alpha=0.5)

# Nombre de frames à conserver dans la trajectoire
trail_length = 50

# Liste pour stocker les trajectoires
trajectories = [([], []) for _ in range(num_objects)]
# print(trajectories)

# Fonction d'initialisation de l'animation
def init():
    points.set_data([], [])
    trails.set_data([], [])
    return points, trails

# Fonction de mise à jour pour chaque frame de l'animation
# def update(frame):
#     x_values = data[:, frame * 2::2]
#     y_values = data[:, frame * 2 + 1::2]

#     points.set_data(x_values, y_values)

#     for i in range(num_objects):
#         # Conserver les 50 dernières valeurs dans la trajectoire
#         trajectories[i] = (trajectories[i][0][-trail_length:] + list(x_values[:, i]),
#                            trajectories[i][1][-trail_length:] + list(y_values[:, i]))

#         trails.set_data(trajectories[i][0], trajectories[i][1])

#     return points, trails

# Nombre total de frames (une frame pour chaque instant de temps)
num_frames = len(data) // 2

# Interval entre les frames en millisecondes
interval = 100

# Créer l'animation
# ani = animation.FuncAnimation(fig, update, frames=num_frames, init_func=init, blit=True, interval=interval)

# Afficher l'animation
# plt.show()
