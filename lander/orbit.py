# Assignment 2: numerical dynamics in three dimensions
import numpy as np
import matplotlib.pyplot as plt


GRAVITY = 6.673E-11 # (m^3/kg/s^2)
MARS_MASS = 6.42E23 # (kg)
MARS_RADIUS = 3386000.0 # (m)
MASS = 100 # (kg)




# Scenario 1: Straight down descent
# Setting up the model as placing the lander on the x-axis, while the centre of Mars lies on the origin
# Plotting function of time


def acceleration(r):
    return - GRAVITY * MARS_MASS * r / np.linalg.norm(r)**3

def simulate(initial_position, initial_velocity, t_max, dt):
    t_array = np.arange(0, t_max, dt)
    N = len(t_array)
    position = np.zeros((N, 3))
    velocity = np.zeros((N, 3))
    position[0], velocity[0] = initial_position, initial_velocity
    position[1] = position[0] + 0.5 * dt**2 * acceleration(position[0]) + dt * velocity[0]
    n_steps = N


    for i in range(1, N-1):
        position[i+1] = 2*position[i] - position[i-1] + dt**2 * acceleration(position[i]) # Verlet implementation of position
        velocity[i] = (position[i+1] - position[i-1]) / (2*dt)
        if np.linalg.norm(position[i+1]) < MARS_RADIUS:
            n_steps = i+1
            break

    velocity[N-1] = (position[N-1] - position[N-2]) / dt

    return t_array[:n_steps], position[:n_steps], velocity[:n_steps]

t_1, position_1, velocity_1 = simulate([MARS_RADIUS + 1E6, 0, 0], [0, 0, 0], 1E3, 1)
altitude = np.linalg.norm(position_1, axis=1) - MARS_RADIUS

# plot the altitude as a function of time
plt.figure(1)
plt.clf()
plt.xlabel('time (s)')
plt.ylabel('altitude (m)')
plt.grid()
plt.plot(t_1, altitude, label='altitude (m)')
plt.legend()
plt.title("Straight down descent")


# Circular orbit
t_2, position_2, velocity_2 = simulate([MARS_RADIUS + 1E6, 0, 0], [0, np.sqrt(GRAVITY * MARS_MASS / (MARS_RADIUS + 1E6)), 0], 1E4, 1)

plt.figure(2)
plt.clf()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.axis('equal')
plt.grid()
plt.plot(position_2[:, 0], position_2[:, 1], label='distance (m)')
plt.legend()
plt.title("Circular orbit")

# Elliptic orbit (initial speed greater than orbital speed, but smaller than escape speed)
t_3, position_3, velocity_3 = simulate([MARS_RADIUS + 1E6, 0, 0], [0, 1E3 + np.sqrt(GRAVITY * MARS_MASS / (MARS_RADIUS + 1E6)), 0], 7 * 1E4, 1)

plt.figure(3)
plt.clf()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.axis('equal')
plt.grid()
plt.plot(position_3[:, 0], position_3[:, 1], label='distance (m)')
plt.legend()
plt.title("Elliptical orbit")


# Hyperbolic escape (initial speed greater than escape speed)
t_4, position_4, velocity_4 = simulate([MARS_RADIUS + 1E6, 0, 0], [0, 2 * np.sqrt(GRAVITY * MARS_MASS / (MARS_RADIUS + 1E6)), 0], 1E4, 1)

plt.figure(4)
plt.clf()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.axis('equal')
plt.grid()
plt.plot(position_4[:, 0], position_4[:, 1], label='distance (m)')
plt.legend()
plt.title("Hyperbolic escape")









plt.show()