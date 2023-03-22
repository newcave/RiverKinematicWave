import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def kinematic_wave_equation(B, L, So, n, Q0, dx, dt, Qin):
    # Compute number of spatial and temporal steps
    nx = int(np.ceil(L/dx))
    nt = len(Qin)

    # Initialize arrays
    h = np.zeros((nt, nx))  # flow depth in ft
    Q = np.zeros((nt, nx))  # flow rate in cfs

    # Initial conditions
    Q[0, :] = Q0
    h[0, :] = Q[0, :]/(B*So)**0.5/n**0.6  # uniform flow

    # Compute the Courant number
    C = dt/(dx/(B*h[0, :])**0.5)

    # Loop over time steps
    for i in range(1, nt):
        # Compute the slope
        S = np.gradient(h[i-1, :])/dx + So

        # Compute the flow rate
        Q[i, :] = Q[i-1, :] - C*B*h[i-1, :]**2*S

        # Compute the flow depth
        h[i, :] = (Q[i, :]/(B*S**0.5*n**0.6))**(3/5)

        # Apply boundary conditions
        h[i, 0] = h[i, 1]
        Q[i, 0] = Q[i, 1]
        h[i, -1] = h[i, -2]
        Q[i, -1] = Q[i, -2]

    # Return results
    return Q

# Constants and parameters
B = 200.0 # channel width in ft
L = 24000.0 # channel length in ft
So = 0.01 # bed slope
n = 0.035 # Manning's roughness factor
Q0 = 2000.0 # initial flow rate in cfs
dx = 3000.0 # spatial step size in ft
dt = 180.0 # temporal step size in sec

# Inflow hydrograph
t = np.arange(0, 7200, dt) # time vector in sec
Qin = np.zeros_like(t) # inflow vector in cfs
Qin[0:20] = 200.0
Qin[20:40] = 400.0
Qin[40:80] = 800.0
Qin[80:160] = 1600.0
Qin[160:320] = 3200.0
Qin[320:640] = 1600.0
Qin[640:720] = 800.0
Qin[720:900] = 400.0

# Run the model
Q = kinematic_wave_equation(B, L, So, n, Q0, dx, dt, Qin)

# Plot the results
X = np.linspace(0, L, Q.shape[1])
T, X = np.meshgrid(t/3600, X)
fig, ax = plt.subplots()
ax.pcolormesh(T, X, Q, cmap='coolwarm')
ax.set_xlabel('Time (hrs.)')
ax.set_ylabel('Distance (meter)')
ax.set_title('Flow rate (CMS)')
plt.colorbar()
st.pyplot(fig)


