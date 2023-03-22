import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def kinematic_wave(B, L, So, n, Q0, dx, dt, inflow):
    nx = int(np.ceil(L/dx))
    nt = len(inflow)
    h = np.zeros((nt, nx))
    Q = np.zeros((nt, nx))
    Q[0, :] = Q0
    h[0, :] = Q[0, :] / (B * So)**0.5 / n**0.6

    C = dt / (dx / (B * h[0, :])**0.5)

    for i in range(1, nt):
        S = np.gradient(h[i-1, :]) / dx + So
        Q[i, :] = Q[i-1, :] - C * B * h[i-1, :]**2 * S
        h[i, :] = (Q[i, :] / (B * S**0.5 * n**0.6))**(3/5)
        h[i, 0] = h[i, 1]
        Q[i, 0] = Q[i, 1]
        h[i, -1] = h[i, -2]
        Q[i, -1] = Q[i, -2]

    return Q, h

B = 100.0
L = 24000.0
So = 0.01
n = 0.025
Q0 = 2000.0
dx = 500.0
dt = 50.0


# Create sliders for inflow values at different time intervals
inflow_200 = st.slider("Inflow for t=0 to 200 (cms)", 0, 1000, 200)
inflow_400 = st.slider("Inflow for t=20 to 400 (cms)", 0, 1000, 400)
inflow_800 = st.slider("Inflow for t=40 to 800 (cms)", 0, 1000, 800)
inflow_1600 = st.slider("Inflow for t=80 to 1600 (cms)", 0, 1000, 1600)
inflow_3200 = st.slider("Inflow for t=160 to 3200 (cms)", 0, 1000, 3200)
inflow_6400 = st.slider("Inflow for t=320 to 6400 (cms)", 0, 1000, 1600)
inflow_7200 = st.slider("Inflow for t=640 to 7200 (cms)", 0, 1000, 800)
inflow_9000 = st.slider("Inflow for t=720 to 9000 (cms)", 0, 1000, 400)

t = np.arange(0, 10000, dt)
inflow = np.zeros_like(t)

inflow = np.zeros_like(t)
inflow[:200] = inflow_200
inflow[200:400] = inflow_400
inflow[400:800] = inflow_800
inflow[800:1600] = inflow_1600
inflow[1600:3200] = inflow_3200
inflow[3200:6400] = inflow_6400
inflow[6400:7200] = inflow_7200
inflow[7200:9000] = inflow_9000


Q, h = kinematic_wave(B, L, So, n, Q0, dx, dt, inflow)

X = np.linspace(0, L, Q.shape[1])
T, X = np.meshgrid(t/3600, X)
fig, axs = plt.subplots(2, 2, figsize=(25, 10))
fig.suptitle('Kinematic Wave Results', fontsize=16)


# Plot the inflow and input values on the same axis (axs[0, 0])
axs[0, 0].plot(t/3600, inflow, label='Inflow')
axs[0, 0].plot(t/3600, inflow, 'ro', label='Input Values')
axs[0, 0].set_xlabel('Time (hours)')
axs[0, 0].set_ylabel('Flow rate (cfs)')
axs[0, 0].legend()

axs[0, 1].plot(X, Q[-1, :])
axs[0, 1].set_xlabel('Distance (meters)')
axs[0, 1].set_ylabel('Flow rate (CMS)')

mesh = axs[1, 0].pcolormesh(T, X, Q[:-1, :-1].T, cmap='coolwarm', shading='flat')
axs[1, 0].set_xlabel('Time (hours)')
axs[1, 0].set_ylabel('Distance (meters)')

# Display the plot in Streamlit
st.pyplot(fig)
