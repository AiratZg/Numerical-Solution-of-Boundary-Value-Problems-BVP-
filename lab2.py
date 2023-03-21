"""
Метод коллокаций
Дифференциальное уравнение: y'' + y' - (4 / x)y = 12x^2 - 2x + 2
Граничные условия: y(1) + y'(1) = 8 && y(2) - y'(2) = 16
Точное решение: y = x^2(x^2 + 1)
Исходный отрезок: [1, 2]
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import diff, simplify, Symbol, factorial, lambdify
t = Symbol('t')

# Исходные данные для задачи
n = 8
alpha_1 = 1
alpha_2 = -1
beta_1 = -1
beta_2 = 1
alpha = 8
beta = 16
a = 1
b = 2
x = (b + a + t * (b - a)) / 2
x_diff = diff(x, t)
p = simplify(1 / (x_diff ** 2))
q = 1 / x_diff
r = simplify(-4 / x)
f = simplify(12 * x**2 - 2 * x + 2)

# Вычисление
alpha_2_new = alpha_2 * 2 / (b - a)
beta_2_new = beta_2 * 2 / (b - a)
d_1_for_z = (beta * alpha_1 - alpha * beta_1) / ((
    alpha_1 + alpha_2_new) * beta_1 + (beta_1 + beta_2_new) * alpha_1)
d_2_for_z = (alpha * (beta_1 + beta_2_new) + beta * (
    alpha_1 + alpha_2_new)) / (
    alpha_1 * (beta_1 + beta_2_new) + beta_1 * (alpha_1 + alpha_2_new))
c_1_for_omega_1 = -(beta_1 * (
    3 * alpha_2_new + alpha_1) + alpha_1 * (
    3 * beta_2_new + beta_1)) / (
        beta_1 * (alpha_1 + alpha_2_new) + alpha_1 * (
                beta_1 + beta_2_new))
d_1_for_omega_1 = ((3 * alpha_2_new + alpha_1) * (
    beta_1 + beta_2_new) - (alpha_1 + alpha_2_new) * (
        3 * beta_2_new + beta_1)) / (alpha_1 * (
                beta_1 + beta_2_new) + beta_1 * (alpha_1 + alpha_2_new))
c_2_for_omega_2 = (beta_1 * (2 * alpha_2_new +
        alpha_1) - alpha_1 * (beta_1 + 2 * beta_2_new)) / (
        alpha_1 * (beta_1 + beta_2_new)
        + beta_1 * (alpha_1 + alpha_2_new))
d_2_for_omega_2 = -((beta_1 + beta_2_new) * (
    2 * alpha_2_new + alpha_1) + (
        beta_1 + 2 * beta_2_new) * (
        alpha_2_new + alpha_1)) / (alpha_1 * (
                beta_1 + beta_2_new) + beta_1 * (alpha_1 + alpha_2_new))
z = d_1_for_z * t + d_2_for_z
L_z = p * diff(z, t, 2) + q * diff(z, t) + r * z
F = simplify(f - L_z)
W = [t ** 3 + c_1_for_omega_1 * t + d_1_for_omega_1,
     t ** 2 + c_2_for_omega_2 * t + d_2_for_omega_2]
grid = [math.cos((2 * i - 1) / (2 * n) * math.pi) for i in range(1, n)]
for i in range(3, n):
    W.append(
        simplify((-1) ** (i - 3) / (2 ** (i - 3) * factorial(i - 3)) * diff(
                (1 - t ** 2) ** (i + 2 - 3), t, i - 3)))
L_i_j = [[float((p * diff(omega, t, 2) + q * diff(omega, t, 1) +
                 r * omega).subs(t, nodes)) for omega in W] for nodes in
         grid]
right_part = [float(F.subs(t, i)) for i in grid]
factors = np.linalg.solve(np.array(L_i_j), np.array(right_part))
y = z
for omega, factor in zip(W, factors):
    y += omega * factor
z = Symbol('z')
analytical = [(x**2 * (x**2+1)).subs(t, (2 * z - b - a) / (
    b - a)).subs(z, (nodes * (b - a) + b + a) / 2) for nodes in grid]
collocation_points = [(nodes * (b - a) + b + a) / 2 for nodes in grid]
xx = np.linspace(a, b, 1000)
yy = y.subs(t, (2 * z - b - a) / (b - a))
approximate_points = [yy.subs(z, nodes) for nodes in collocation_points]
analytical_pow_2 = [item ** 2 for item in analytical]
Q_ocm = [(first - second) ** 2 for first, second in zip(
    analytical, approximate_points)]
numerator = sum(Q_ocm)
denominator = sum(analytical_pow_2)


def builder():
    """Построение графика"""
    approximate = lambdify(z, yy)(xx)
    plt.plot(xx, np.transpose(approximate), color='black')
    plt.plot(collocation_points, analytical, 'o', color='red')
    plt.show()


def error():
    """Вычисление ошибки"""
    print(f"Error: {numerator / denominator * 100}%")


if __name__ == '__main__':
    builder()
    print(f"Приближенное решение: {factors};")
    error()
