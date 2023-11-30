import pandas as pd
import csv
import sys

def heal(in_filename,out_filename,nb_simulations):

    #nb_simulations = 3
    #file_name = "Cured_Triangle.csv"

    # Charger le fichier CSV
    data = pd.read_csv(in_filename, delimiter=";", header=None)

    # Nombre de lignes du fichier final
    final_lines = int(data.shape[0]/nb_simulations)

    # Vérifier la forme initiale du DataFrame
    print("Forme initiale du DataFrame :", data.shape)

    # Créer le nouveau fichier CSV
    new_file = pd.DataFrame().to_csv(out_filename, index=False, header=False)

    with open(out_filename , 'w', encoding='UTF8', newline='') as file:

        writer = csv.writer(file, delimiter=';')
        for i in range(final_lines):
            
            recipe = []
            for sim in range(nb_simulations):
                recipe.append(data.iloc[i+sim*final_lines])
                
            writer.writerow(pd.concat(recipe))

    new_file = pd.read_csv(out_filename, delimiter=";", header=None)

    # Vérifier la forme du nouveau DataFrame
    print("Forme du nouveau DataFrame :", new_file.shape)

nb_simulations = int(sys.argv[1])

heal('Triangle.csv','Cured_Triangle.csv',nb_simulations)
heal('TriangleV.csv','Cured_Triangle_v.csv',nb_simulations)