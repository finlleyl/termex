import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Константы системы
R = 5.0  # радиус большой шестерни
r = 1.5  # радиус маленькой шестерни
g = 9.81  # ускорение свободного падения
m1, m2, m3 = 10.0, 5.0, 2.0  # массы
c = 50.0  # жесткость пружины
M1, M2 = 0.0, 0.0  # моменты

# Переменные
phi, phi_dot, phi_ddot = sp.symbols('phi phi_dot phi_ddot')
theta, theta_dot, theta_ddot = sp.symbols('theta theta_dot theta_ddot')

# Лагранжиан
T_rot1 = (1 / 2) * (m1 + m2) * (R**2) * theta_dot**2
T_rot2 = (1 / 2) * m2 * (R - r)**2 * phi_dot**2
T_crank = (1 / 6) * m3 * (R - r)**2 * phi_dot**2
T = T_rot1 + T_rot2 + T_crank

V_potential = m2 * g * (R - r) * sp.sin(phi) + (1 / 2) * c * phi**2
V = V_potential

L = T - V

# Уравнения Лагранжа
L_theta = sp.diff(sp.diff(L, theta_dot), 't') - sp.diff(L, theta)
L_phi = sp.diff(sp.diff(L, phi_dot), 't') - sp.diff(L, phi)

# Упрощение уравнений
eq_theta = sp.simplify(L_theta).subs('t', 0)  # Убираем временные производные
eq_phi = sp.simplify(L_phi).subs('t', 0)

# Решение системы уравнений
a11 = eq_phi.coeff(phi_ddot, 1)
a12 = eq_phi.coeff(theta_ddot, 1)
a21 = eq_theta.coeff(phi_ddot, 1)
a22 = eq_theta.coeff(theta_ddot, 1)
b1 = -eq_phi.subs([(phi_ddot, 0), (theta_ddot, 0)])
b2 = -eq_theta.subs([(phi_ddot, 0), (theta_ddot, 0)])

# Числовые функции
phi_ddot_func = sp.lambdify((phi, phi_dot, theta, theta_dot), a11 * a22 - a12 * a21)
theta_ddot_func = sp.lambdify((phi, phi_dot, theta, theta_dot), b2 - b1)

# Численное решение системы
def equations(t, y):
    phi, phi_dot, theta, theta_dot = y
    phi_ddot_val = phi_ddot_func(phi, phi_dot, theta, theta_dot)
    theta_ddot_val = theta_ddot_func(phi, phi_dot, theta, theta_dot)
    return [phi_dot, phi_ddot_val, theta_dot, theta_ddot_val]

# Начальные условия
y0 = [0, 0, 0, 0]  # phi, phi_dot, theta, theta_dot
t_span = (0, 10)
t_eval = np.linspace(*t_span, 1000)

solution = solve_ivp(equations, t_span, y0, t_eval=t_eval)

# Построение графиков
phi_values = solution.y[0]
theta_values = solution.y[2]

plt.figure(figsize=(10, 6))
plt.plot(t_eval, phi_values, label=r'$\phi(t)$')
plt.plot(t_eval, theta_values, label=r'$\theta(t)$')
plt.xlabel('Время (с)')
plt.ylabel('Углы (рад)')
plt.title('Движение системы')
plt.legend()
plt.grid()
plt.show()
