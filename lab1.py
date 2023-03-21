"""Метод прогонки."""
"""
Дифференциальное уравнение: y'' + y'(x-1) - 3y = 6(x-1) + (x-3) * e^x
Граничные условия: y(0) + y'(0) = 4 && y(1) + y'(1) = 2e
Точное решение: y = x^2(x^2 + 1)
Исходный отрезок: [0, 1]
"""


import math
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from numpy import arange

flag_ten = 1
while (True):
    if flag_ten:
        print(
            "Краевая задача вида: \n -p(x) * y'' + q(x) * y' + r(x) * y = f(x), a<=x<=b, "
            "\n alpha_1 * y(a) - alpha_2 * y'(a) = alpha, "
            "\n beta_1 * y(b) + beta_2 * y'(b) = beta. Численное решение представляет порядок точности O(h).")
        ALPHA = 4.0
        BETA = 2.0 * math.exp(1)
        N = 8
        b = 1.0
        a = 0.0
        h = (b - a) / N
        p = -1.0
        q = 1.0
        r = -3.0
        ALPHA_1 = 1.0
        ALPHA_2 = -1.0
        BETA_1 = 1.0
        BETA_2 = 1.0
        exact_solution_temp = np.linspace(a - h / 2, b + h / 2, 100)
        EXACT_SOLUTION = [(x - 1) ** 3 + math.exp(x) for x in exact_solution_temp]
        X_I_extended = list(arange(a - h / 2, b + h / 2 + h, h))
        P_I = [p for x in X_I_extended]
        Q_I = [x - 1 for x in X_I_extended]
        R_I = [r for x in X_I_extended]
        F_I = [6 * (x - 1) + (x - 3) * math.exp(x) for x in X_I_extended]
        A = [0.0]
        C = [ALPHA_1 / 2 - ALPHA_2 / h]
        B = [-(ALPHA_1 / 2 + ALPHA_2 / h)]
        s = [C[0] / B[0]]
        G = [ALPHA]
        t = [-G[0] / B[0]]
        y_exact = [(x - 1) ** 3 + math.exp(x) for x in X_I_extended]
        flag_ten = True
        break
    else:
        break


def print_table(field_names, row):
    """Распечатаем итоговую таблицу"""
    table = PrettyTable()
    table.add_column("i", [i for i in range(0, len(row[0]))])
    for i, field in enumerate(field_names):
        table.add_column(field, [round(item, 5) for item in row[i]])
    print(table)


def sum_square(approximate, exact):
    """Посчитаем остаточную сумму квадратов"""
    return [(vi - wi) ** 2 for vi, wi in zip(approximate, exact)]


def error(y_i_approximate):
    """Погрешность по среднеквадратичным отклонениям"""
    print("Оценка погрешности")
    remainder = sum_square(y_exact, y_i_approximate)
    sum_square_y = [i * i for i in y_exact]
    print_table(["x_i", "y_i", "y_i_[approximation]", "(y_i - y_i_[approximation])^2", "y_i^2"],
                [X_I_extended, y_exact, y_i_approximate, remainder, sum_square_y])
    remainder = sum(remainder)
    sum_square_y = sum(sum_square_y)
    print(f"Погрешность {round(remainder / sum_square_y * 100, 4)} %")


def build_graph(graph_name, y_i_approximate):
    """Строим графики"""
    exact_solution_array = np.linspace(a - h / 2, b + h / 2, 100)
    if flag_ten:
        schedule_exact = EXACT_SOLUTION
    else:
        schedule_exact = [eval(EXACT_SOLUTION) for _ in exact_solution_array]
    plt.title(graph_name)
    plt.plot(exact_solution_array, schedule_exact)
    plt.scatter(X_I_extended, y_i_approximate)
    plt.show()


def second_order_approximation():
    """Аппроксимация второго порядка, метод прогонки"""
    for i in range(1, N + 1):
        A.append(-P_I[i] / h ** 2 - Q_I[i] / (2 * h))
        C.append(-P_I[i] / h ** 2 + Q_I[i] / (2 * h))
        B.append(-2 * P_I[i] / h ** 2 - R_I[i])
    A.append(BETA_1 / 2 - BETA_2 / h)
    C.append(0.0)
    B.append(-(BETA_1 / 2 + BETA_2 / h))
    for i in range(1, N + 1):
        s.append(C[i] / (B[i] - A[i] * s[i - 1]))
    s.append(0.0)
    for f_i in F_I[1:N + 1]:
        G.append(f_i)
    G.append(BETA)
    for i in range(1, N + 2):
        t.append((A[i] * t[i - 1] - G[i]) / (B[i] - A[i] * s[i - 1]))
    temp = 0
    y_approximate = [t[-1]]
    for s_i, t_i in zip(s[::-1][1:], t[::-1][1:]):
        y_approximate.append(s_i * y_approximate[temp] + t_i)
        temp = temp + 1
    y_approximate = y_approximate[::-1]
    print("Расчет оформим в виде следующей таблицы")
    print_table(["x_i", "A_i", "B_i", "C_i", "G_i", "s_i", "t_i", "y_i"],
                [X_I_extended, A, B, C, G, s, t, y_approximate])
    error(y_approximate)
    build_graph("Аппроксимация O(h^2)", y_approximate)


if __name__ == '__main__':
    second_order_approximation()
