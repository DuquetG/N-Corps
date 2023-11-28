import plotly.graph_objects as go
import numpy as np

# Create data for two objects
t = np.linspace(0, 10, 100)
x1 = np.sin(t)
y1 = np.cos(t)
x2 = np.cos(t)
y2 = np.sin(t)

# Create figure
fig = go.Figure()

# Add scatter trace for the first marker
scatter_trace1 = go.Scatter(
    x=[x1[0]],
    y=[y1[0]],
    mode='markers',
    marker=dict(size=10, color='red'),
    name='Marker 1'
)
fig.add_trace(scatter_trace1)

# Add line trace for the connecting line of the first object
line_trace1 = go.Scatter(
    x=[x1[0]],
    y=[y1[0]],
    mode='lines',
    line=dict(color='blue', width=2),
    name='Line 1'
)
fig.add_trace(line_trace1)

# Add scatter trace for the second marker
scatter_trace2 = go.Scatter(
    x=[x2[0]],
    y=[y2[0]],
    mode='markers',
    marker=dict(size=10, color='green'),
    name='Marker 2'
)
fig.add_trace(scatter_trace2)

# Add line trace for the connecting line of the second object
line_trace2 = go.Scatter(
    x=[x2[0]],
    y=[y2[0]],
    mode='lines',
    line=dict(color='orange', width=2),
    name='Line 2'
)
fig.add_trace(line_trace2)

# Set up animation
frames = [go.Frame(data=[
    go.Scatter(x=x1[:i+1], y=y1[:i+1], mode='lines', line=dict(color='blue', width=2)),
    go.Scatter(x=[x1[i]], y=[y1[i]], mode='markers', marker=dict(size=10, color='red')),
    go.Scatter(x=x2[:i+1], y=y2[:i+1], mode='lines', line=dict(color='orange', width=2)),
    go.Scatter(x=[x2[i]], y=[y2[i]], mode='markers', marker=dict(size=10, color='green'))
]) for i in range(1, len(t))]

# print(x1, x1[10+1], [x1[10]])

fig.frames = frames

# Update layout
fig.update_layout(
    updatemenus=[
        dict(
            type='buttons',
            showactive=False,
            buttons=[dict(label='Play',
                          method='animate',
                          args=[None, dict(frame=dict(duration=100, redraw=True), fromcurrent=True)])]
        )
    ]
)

# Update layout for the initial state
fig.update_layout(
    xaxis=dict(range=[-2, 2]),
    yaxis=dict(range=[-2, 2]),
    title='Animated Markers and Lines for Two Objects'
)

fig.show()
