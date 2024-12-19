import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

R = 5.0  # радиус большой шестерни
r = 1.5  # радиус маленькой шестерни
omega = 0.5  # угловая скорость

iterations = 200
T = np.linspace(0, 10, iterations)


def theta(t):
    return omega * t


def phi(t):
    return 2.0 * omega * t


# Позиции точки A:
# A движется по окружности радиуса (R - r) вокруг O
# X_A(t) = (R - r)*cos(phi(t))
# Y_A(t) = (R - r)*sin(phi(t))

XA = (R - r) * np.cos(phi(T))
YA = (R - r) * np.sin(phi(T))


def draw_circle(center_x, center_y, radius=1.0, resolution=100):
    angles = np.linspace(0, 2 * math.pi, resolution)
    x = center_x + radius * np.cos(angles)
    y = center_y + radius * np.sin(angles)
    return x, y


fig, ax = plt.subplots()
ax.set_aspect('equal', adjustable='box')
ax.set_xlim(-R - 2, R + 2)
ax.set_ylim(-R - 2, R + 2)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Анимация системы шестерен')

# (Шестерня 1)
gear1_x, gear1_y = draw_circle(0, 0, radius=R)
gear1_line, = ax.plot(gear1_x, gear1_y, 'b')

# (Шестерня 2)
gear2_x, gear2_y = draw_circle(XA[0], YA[0], radius=r)
gear2_line, = ax.plot(gear2_x, gear2_y, 'r')

# кривошип OA
crank_line, = ax.plot([0, XA[0]], [0, YA[0]], 'k')

point_O, = ax.plot(0, 0, 'ko')
point_A, = ax.plot(XA[0], YA[0], 'ro')


def update(i):
    gear2_x, gear2_y = draw_circle(XA[i], YA[i], radius=r)
    gear2_line.set_data(gear2_x, gear2_y)

    crank_line.set_data([0, XA[i]], [0, YA[i]])

    point_A.set_data([XA[i]], [YA[i]])

    return gear2_line, crank_line, point_A



anim = FuncAnimation(fig, update, frames=iterations, interval=50, blit=True)

plt.show()
