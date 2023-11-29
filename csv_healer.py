import pandas as pd
import csv

nb_simulations = 100
file_name = "Cured_Triangle.csv"

# Charger le fichier CSV
data = pd.read_csv('Triangle.csv', delimiter=";", header=None)

# Nombre de lignes du fichier final
final_lines = int(data.shape[0]/nb_simulations)

# Vérifier la forme initiale du DataFrame
print("Forme initiale du DataFrame :", data.shape)

print(range(nb_simulations))

# Créer le nouveau fichier CSV
new_file = pd.DataFrame().to_csv(file_name, index=False, header=False)

with open(file_name , 'w', encoding='UTF8', newline='') as file:

    writer = csv.writer(file, delimiter=';')
    for i in range(final_lines):
        
        recipe = []
        for sim in range(nb_simulations):
            recipe.append(data.iloc[i+sim*final_lines])
            
        writer.writerow(pd.concat(recipe))

new_file = pd.read_csv(file_name, delimiter=";", header=None)

# Vérifier la forme du nouveau DataFrame
print("Forme du nouveau DataFrame :", new_file.shape)
