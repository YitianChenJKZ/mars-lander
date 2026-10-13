# uncomment the next line if running in a notebook
# %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

# mass, spring constant, initial position and velocity
m = 1
k = 1
x = 0
v = 1

# simulation time, timestep and time
t_max = 300
dt = 1 # Critical value of dt for Verlet algorithm to become unstable is dt = 1
t_array = np.arange(0, t_max, dt)

# initialise empty lists to record trajectories
x_list_euler = []
v_list_euler = []
x_list_verlet = []
v_list_verlet = []

# Euler integration
for t in t_array:

    # append current state to trajectories
    x_list_euler.append(x)
    v_list_euler.append(v)

    # calculate new position and velocity
    a = -k * x / m
    x = x + dt * v
    v = v + dt * a


# mass, spring constant, initial position and velocity
m = 1
k = 1
x_current = 0
v = 1
x_previous = x_current - v*dt



# Verlet integration
for t in t_array:

    # append current state to trajectories
    x_list_verlet.append(x_current)
    v_list_verlet.append(v)

    # calculate new position and velocity
    a = -k * x_current / m
    x_new = a * dt**2 + 2 * x_current - x_previous
    v = (x_new - x_previous) / (2*dt)
    x_previous = x_current
    x_current = x_new


# convert trajectory lists into arrays, so they can be sliced (useful for Assignment 2)
x_array_euler = np.array(x_list_euler)
v_array_euler = np.array(v_list_euler)
x_array_verlet = np.array(x_list_verlet)
v_array_verlet = np.array(v_list_verlet)











# plot the position-time graph
plt.figure(1)
plt.clf()
plt.xlabel('time (s)')
plt.grid()
plt.plot(t_array, x_array_euler, label='x (m)')
plt.plot(t_array, v_array_euler, label='v (m/s)')
plt.legend()

plt.figure(2)
plt.clf()
plt.xlabel('time (s)')
plt.grid()
plt.plot(t_array, x_array_verlet, label='x (m)')
plt.plot(t_array, v_array_verlet, label='v (m/s)')
plt.legend()


plt.show()
