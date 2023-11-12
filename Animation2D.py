import pandas as pd
import plotly.graph_objects as go

# Charger le fichier CSV avec pandas
data = pd.read_csv("bodies_movement2D.csv", delimiter=";")

# Créer une figure Plotly
fig = go.Figure()

# Ajouter une trace pour chaque objet
for i in range(0, len(data.columns), 2):
    trace = go.Scatter(x=data.iloc[:, i], y=data.iloc[:, i + 1],
                      mode='markers', name=f'Objet {i//2 + 1}')
    fig.add_trace(trace)

# Configurer les paramètres de la mise en page
fig.update_layout(title_text="Positions des objets au fil du temps",
                  xaxis_title="Position X", yaxis_title="Position Y",
                  showlegend=True)

# Ajouter les frames pour l'animation
frames = [go.Frame(data=[go.Scatter(x=data.iloc[:frame + 1, i], y=data.iloc[:frame + 1, i + 1],
                                   mode='markers', name=f'Objet {i//2 + 1}') for i in range(0, len(data.columns), 2)]) for frame in range(1, len(data))]

# Ajouter les frames à la figure
fig.frames = frames

# Configurer les paramètres de l'animation
animation_settings = dict(frame=dict(duration=1, redraw=True), fromcurrent=True)
fig.update_layout(updatemenus=[dict(type='buttons', showactive=False,
                                    buttons=[dict(label='Play',
                                                  method='animate',
                                                  args=[None, animation_settings])])])

# Afficher le graphique interactif
fig.show()
