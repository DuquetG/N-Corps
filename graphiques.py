import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

# Charger le fichier CSV avec pandas
#data = pd.read_csv("bodies_movement2D.csv", delimiter=";", skiprows=lambda x: x % 1000 != 0)
#data=data.to_numpy()
#Graphique des positions dans le temps, avec code couleur pour chaque masse
def plot_positions():
    data = pd.read_csv("bodies_movement2D.csv", delimiter=";", skiprows=lambda x: x % 1000 != 0)
    data=data.to_numpy()
    for i in range(0,math.floor(data.shape[1]/2)):
        plt.plot(data[:,2*i],data[:,2*i+1])
    plt.title('Positions au fil du temps')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

#Graphique de l'énergie en fonction du temps, énergie cinétique, potentielle puis la somme des deux
def plot_energy():
    data = pd.read_csv("energy.csv", delimiter=";", skiprows=lambda x: x % 1000 != 0)
    data=data.to_numpy()
    for i in range(0,math.floor(data.shape[1]-1)):
        plt.plot(data[:,-1],data[:,i])
    plt.title('Énergies au fil du temps')
    plt.xlabel("t")
    plt.ylabel("E")
    plt.show()



#Graphique de l'énergie potentielle et cinétique moyenne en fonction du temps
def plot_mean_energies():
    data = pd.read_csv("viriel.csv", delimiter=";", skiprows=lambda x: x % 1000 != 0)
    data=data.to_numpy()
    for i in range(0,math.floor(data.shape[1]-1)):
        plt.plot(data[:,-1],data[:,i])
    plt.title('Énergies moyennes au fil du temps')
    plt.xlabel("t")
    plt.ylabel("<E>")
    plt.show()

#Graphique de la distance entre 2 masses en fonction du temps (on devrait avoir une courbe exponentielle)

#Graphique de l'exposant de Lyapunov
def plot_Lyapunov():
    data = pd.read_csv("Lyapunov.csv", delimiter=";",skiprows=lambda x: x % 10 != 0)
    data=data.to_numpy()
    legende = []
    for i in range(1,math.floor(data.shape[1])):
        legende.append(f"L{i}")
        plt.plot(data[:,0],data[:,i])
    plt.axhline(y=0, color='black', linestyle='--', label='y=0')
    plt.plot()
    plt.title('Exposant de Lyapunov au fil du temps')
    plt.xlabel("t")
    plt.ylabel("Lambda")
    plt.legend((legende))

#Graphique temps de Lyapunov
def plot_Lyapunov_time():
    data = pd.read_csv("Lyapunov.csv", delimiter=";",skiprows=lambda x: x % 10 != 0)
    data=data.to_numpy()
    legende = []
    for i in range(1,math.floor(data.shape[1])):
        legende.append(f"T{i}")
        plt.plot(data[:,0],1/data[:,i])
    #plt.axhline(y=0, color='black', linestyle='--', label='y=0')
    plt.plot()
    plt.title('Temps de Lyapunov au fil du temps')
    plt.xlabel("t")
    plt.ylabel("Tau")
    plt.legend((legende))
#Exécution des graphiques

plot_Lyapunov_time()
plt.show()