import streamlit as st
import numpy as np
import pandas as pd
import math 
import matplotlib.pyplot as plt

def kinematic_wave(B, L, So, n, Q0, dx, dt, inflow):
    nx = math.ceil(L/dx)
    nt = len(inflow)

    h = [[0] * nx for i in range(nt)]
    Q = [[0] * nx for i in range(nt)]
    S = [0] * nx

    # initialize Q and h at t = 0
    for i in range(nx):
        Q[0][i] = Q0
        h[0][i] = Q0 / (B * math.pow(So, 0.5) * math.pow(n, 0.6))
        st.write(i,'에서의 h[0][i]값=', h[0][i])

    # calculate C at t = 0
    C = dt / (dx / (B * math.pow(h[0][0], 0.5)))
    st.write(f't = 0, x = 0, C = {C:.4f}')
    st.write('원래의 C값')

    # calculate S at t = 0
    S[0] = (h[0][1] - h[0][0]) / dx + So
    S[nx-1] = (h[0][nx-1] - h[0][nx-2]) / dx + So
    for i in range(1, nx-1):
        S[i] = (h[0][i+1] - h[0][i-1]) / (2 * dx) + So

    # solve the kinematic wave equation for all t > 0
    for i in range(1, nt):

        # calculate Q at t = i
        for j in range(nx):
            Q[i][j] = Q[i-1][j] - C * B * math.pow(h[i-1][j], 2/3) * math.pow(S[j], 1/2) * dt / dx

        # calculate h at t = i
        for j in range(nx):
            h[i][j] = Q[i][j] / (B * math.pow(S[j], 1/2) * math.pow(n, 1/6))

        # enforce boundary conditions
        h[i][0] = h[i][1]
        Q[i][0] = Q[i][1]
        h[i][nx-1] = h[i][nx-2]
        Q[i][nx-1] = Q[i][nx-2]

        # update C and S for the next iteration
        st.write('calc started...')
        st.write(i, h[i][0])
        if h[i][0] > 0:
            C = dt / (dx / (B * math.pow(h[i][0], 0.5)))
        
        st.write(f't = {i*dt:.1f}, x = 0, C = {C:.4f}')
        st.write('updated C value')
        S[0] = (h[i][1] - h[i][0]) / dx + So
        S[nx-1] = (h[i][nx-1] - h[i][nx-2]) / dx + So
        for j in range(1, nx-1):
            S[j] = (h[i][j+1] - h[i][j-1]) / (2 * dx) + So

    # return the final flow rate at the downstream end of the channel
    return Q[nt-1][nx-1]

st.write("| Parameter | 주요 기본변수  |")
st.write("| ----------------------------------------|")


B = st.slider('Channel bottom width (meters)', min_value=1.0, max_value=1000.0, value=100.0)
L = st.slider('Channel length (meters)', min_value=1000.0, max_value=100000.0, value=10000.0, step=1000.0)
So = st.slider('Channel slope', min_value=0.001, max_value=0.1, value=0.05, step=0.01)
n = st.slider('Manning roughness coefficient', min_value=0.01, max_value=0.05, value=0.025, step=0.001)
Q0 = st.slider('Flow rate at upstream boundary (CMS)', min_value=1.0, max_value=500.0, value=200.0)
dx = st.slider('Distance step size (meters)', min_value=10.0, max_value=1000.0, value=500.0)
dt = st.slider('Time step size (seconds)', min_value=60, max_value=3600, value=180)

#B = 100.0
#L = 10000.0
#So = 0.05
#n = 0.025
#Q0 = 200.0
#dx = 500.0
#dt = 180


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
#time_intervals = [4 * 3600, 8 * 3600, 12 * 3600, 16 * 3600, 20 * 3600, 24 * 3600, 28 * 3600, 32 * 3600]
#inflow_values = [inflow_0, inflow_4, inflow_8, inflow_12, inflow_16, inflow_20, inflow_24, inflow_28]

#for i in range(len(time_intervals) - 1):
#    start_time = time_intervals[i] // dt
#    end_time = time_intervals[i + 1] // dt
#    inflow[start_time:end_time] = inflow_values[i]
    
### 인터벌 고려해서 입력 받음 끝

inflow[:4 * 3600 //dt]= inflow_0
inflow[4 * 3600 //dt:8 * 3600 //dt] = inflow_4
inflow[8 * 3600 //dt:12 * 3600 //dt] = inflow_8
inflow[12 * 3600 //dt:16 * 3600 //dt] = inflow_12
inflow[16 * 3600 //dt:20 * 3600 //dt] = inflow_16
inflow[20 * 3600 //dt:24 * 3600 //dt] = inflow_20
inflow[24 * 3600 //dt:28 * 3600 //dt] = inflow_24
inflow[28 * 3600 //dt:32 * 3600 //dt] = inflow_28

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

# 초기 조건에서 강 유량 시각화
axs[0, 1].plot(X, Q[0, :])
axs[0, 1].set_xlabel('Distance (meters)')
axs[0, 1].set_ylabel('Flow rate (CMS)')

axs[1, 0].plot(t, h, label='water level')
axs[1, 0].set_xlabel('Time (seconds)')
axs[1, 0].set_ylabel('Flow rate (CMS)')
axs[1, 0].legend()

t_list = [4*3600, 8*3600, 10*3600] # 4시간, 12시간, 20시간
for i, t_idx in enumerate(t_list):
    axs[1, 1].plot(X, Q[t_idx // dt, :], label=f'{t_idx // 3600} 시간')  # 선택한 시간에 해당하는 인덱스를 사용
    axs[1, 1].set_xlabel('Distance (meters)')
    axs[1, 1].set_ylabel('Flow rate (CMS)')
    axs[1, 1].legend(['After 4hrs', '8hrs', '10hrs'])




#im = axs[1, 0].contourf(X, t, h.T, cmap='YlGnBu')
#axs[1, 0].set_xlabel('Distance (meters)')
#axs[1, 0].set_ylabel('Time (hours)')
#cbar = fig.colorbar(im, ax=axs[0, 0])
#cbar.set_label('Water depth (meters)')

plt.tight_layout()
st.pyplot(fig)


# Save the results to a CSV file for a given time interval
time_interval = 2 * 3600 # 2 hours
result_df = pd.DataFrame({'Distance (meters)': X[0, :], 'Flow rate (CMS)': Q[(time_interval // dt), :]})
result_df.to_csv(f'kinematic_wave_result_{time_interval // 3600}hrs.csv', index=False)

# Plot the results at a given time interval
axs[1, 1].plot(X, Q[time_interval // dt, :], label=f'{time_interval // 3600} hours')
axs[1, 1].set_xlabel('Distance (meters)')
axs[1, 1].set_ylabel('Flow rate (CMS)')
axs[1, 1].legend()

# Save the figure
plt.tight_layout()
st.pyplot(fig)

# Save the CSV file
st.markdown(f"### Result for {time_interval // 3600} hours")
st.write(result_df)
st.markdown(filedownload(result_df, f'kinematic_wave_result_{time_interval // 3600}hrs.csv'), unsafe_allow_html=True)

# Function to create a download link for a given dataframe
def filedownload(df, filename):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download CSV File</a>'
    return href
