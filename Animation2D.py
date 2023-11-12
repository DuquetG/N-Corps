import pandas as pd
import plotly.graph_objects as go

# Charger le fichier CSV avec pandas
data = pd.read_csv("bodies_movement2D.csv", delimiter=";", skiprows=lambda x: x % 1000 != 0)

# Créer une liste de couleurs pour chaque objet
colors = ["lemonchiffon", "silver", "goldenrod", "olivedrab", "indianred", "darksalmon", "sandybrown", "lightsteelblue", "blueviolet"]

# Créer une figure Plotly
fig = go.Figure()

# Ajouter une trace pour chaque objet
for i in range(0, len(data.columns), 2):
    trace = go.Scatter(x=data.iloc[:, i], y=data.iloc[:, i + 1],
                      mode='markers', name=f'Objet {i//2 + 1}',
                      marker=dict(color=colors[i//2]))
    fig.add_trace(trace)

# Configurer les paramètres de la mise en page
fig.update_layout(title_text="Positions des masses au fil du temps",
                  xaxis_title="Position X", yaxis_title="Position Y",
                  showlegend=True,
                  plot_bgcolor='black',
                  paper_bgcolor='black',
                  template="plotly_dark")
#                   xaxis=dict(range=[min(data.iloc[:, ::2].values.flatten()), max(data.iloc[:, ::2].values.flatten())]),
#                   yaxis=dict(range=[min(data.iloc[:, 1::2].values.flatten()), max(data.iloc[:, 1::2].values.flatten())])
# )

# Ajouter les frames pour l'animation
frames = [go.Frame(data=[go.Scatter(x=data.iloc[:frame + 1, i], y=data.iloc[:frame + 1, i + 1],
                                   mode='markers', name=f'Masse {i//2 + 1}') for i in range(0, len(data.columns), 2)]) for frame in range(1, len(data))]

# Ajouter les frames à la figure
fig.frames = frames

# Configurer les paramètres de l'animation
animation_settings = dict(frame=dict(duration=0.1, redraw=True), fromcurrent=True)
fig.update_layout(updatemenus=[dict(type='buttons', showactive=False,
                                    buttons=[dict(label='Play',
                                                  method='animate',
                                                  args=[None, animation_settings])])])

# Afficher le graphique interactif
fig.show()
