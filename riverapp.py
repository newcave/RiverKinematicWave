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
L = 2500.0
So = 0.05
n = 0.025
Q0 = 2000.0
dx = 500.0
dt = 30.0


# Create sliders for inflow values at different time intervals
inflow_range = st.slider("총분석시간(hrs)", 0, 48, 12)
inflow_0 = st.slider("Inflow for t=0 to 4 hrs. ", 0, 5000, 200)
inflow_4 = st.slider("Inflow for t=4 to 8 hrs.", 0, 5000, 400)
inflow_8 = st.slider("Inflow for t=8 to 12 hrs.", 0, 5000, 800)
inflow_12 = st.slider("Inflow for t=12 to 16 hrs.", 0, 5000, 1600)
inflow_16 = st.slider("Inflow for t=16 to 20 hrs.", 0, 5000, 3200)
inflow_20 = st.slider("Inflow for t=20 to 24 hrs.", 0, 5000, 1600)
inflow_24 = st.slider("Inflow for t=24 to 28 hrs.", 0, 5000, 800)
inflow_28 = st.slider("Inflow for t=28 to 32 hrs.", 0, 5000, 400)

t = np.arange(0, inflow_range*3600, dt)

inflow = np.zeros_like(t)

### 인터벌 대로 입력 받음
time_intervals = [4 * 3600, 8 * 3600, 12 * 3600, 16 * 3600, 20 * 3600, 24 * 3600, 28 * 3600, 32 * 3600]
inflow_values = [inflow_0, inflow_4, inflow_8, inflow_12, inflow_16, inflow_20, inflow_24, inflow_28]

for i in range(len(time_intervals) - 1):
    start_time = time_intervals[i] // dt
    end_time = time_intervals[i + 1] // dt
    inflow[start_time:end_time] = inflow_values[i]

    
### 인터벌 고려해서 입력 받음 끝

#inflow[:4 * 3600] = inflow_0
#inflow[4 * 3600:8 * 3600] = inflow_4
#inflow[8 * 3600:12 * 3600] = inflow_8
#inflow[12 * 3600:16 * 3600] = inflow_12
#inflow[16 * 3600:20 * 3600] = inflow_16
#inflow[20 * 3600:24 * 3600] = inflow_20
#inflow[24 * 3600:28 * 3600] = inflow_24
#inflow[28 * 3600:32 * 3600] = inflow_28

Q, h = kinematic_wave(B, L, So, n, Q0, dx, dt, inflow)

X = np.linspace(0, L, Q.shape[1])
T, X = np.meshgrid(t, X)
fig, axs = plt.subplots(2, 2, figsize=(25, 10))
fig.suptitle('Kinematic Wave Results', fontsize=16)


# Plot the inflow and input values on the same axis (axs[0, 0])
axs[0, 0].plot(t, inflow, label='Inflow')
axs[0, 0].set_xlabel('Time (seconds)')
axs[0, 0].set_ylabel('Flow rate (CMS)')
axs[0, 0].legend()

axs[0, 1].plot(t*dt, inflow, 'ro', label='Input Values')
axs[0, 1].set_xlabel('Time (seconds)')
axs[0, 1].set_ylabel('Flow rate (CMS)')
axs[0, 1].legend()


#axs[0, 1].plot(X, Q[-1, :])
#axs[0, 1].set_xlabel('Distance (meters)')
#axs[0, 1].set_ylabel('Flow rate (CMS)')

mesh = axs[1, 0].pcolormesh(T, X, Q[:-1, :-1].T, cmap='coolwarm', shading='flat')
axs[1, 0].set_xlabel('Time (hours)')
axs[1, 0].set_ylabel('Distance (meters)')

# Display the plot in Streamlit
st.pyplot(fig)
