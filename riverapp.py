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

B = 200.0
L = 24000.0
So = 0.01
n = 0.035
Q0 = 2000.0
dx = 3000.0
dt = 180.0

t = np.arange(0, 7200, dt)
inflow = np.zeros_like(t)
inflow[:20] = 200.0
inflow[20:40] = 400.0
inflow[40:80] = 800.0
inflow[80:160] = 1600.0
inflow[160:320] = 3200.0
inflow[320:640] = 1600.0
inflow[640:720] = 800.0
inflow[720:900] = 400.0

Q, h = kinematic_wave(B, L, So, n, Q0, dx, dt, inflow)

X = np.linspace(0, L, Q.shape[1])
T, X = np.meshgrid(t/3600, X)
fig, axs = plt.subplots(2, 2, figsize=(25, 10))
fig.suptitle('Kinematic Wave Results', fontsize=16)

axs[0, 0].plot(t/3600, inflow)
axs[0, 0].set_xlabel('Time (hours)')
axs[0, 0].set_ylabel('Flow rate (cfs)')

axs[0, 1].plot(X, Q[-1, :])
axs[0, 1].set_xlabel('Distance (meters)')
axs[0, 1].set_ylabel('Flow rate (CMS)')

mesh = axs[1, 0].pcolormesh(T, X, Q[:-1, :-1].T, cmap='coolwarm', shading='flat')
axs[1, 0].set_xlabel('Time (hours)')
axs[1, 0].set_ylabel('Distance (meters)')

# Display the plot in Streamlit
st.pyplot(fig)
