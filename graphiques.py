import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Charger le fichier CSV avec pandas
data = pd.read_csv("bodies_movement2D.csv", delimiter=";", skiprows=lambda x: x % 1000 != 0)
