import math
import matplotlib.pyplot as plt

# Задаємо функції
def f_x(t, x, y):
    return 998*x + 1998*y

def f_y(t, x, y):
    return 999*x - 1999*y

# Метод Adams-Bashforth 4-го порядку
def adams_bashforth_4(t_n, x_n, y_n, h):
    t_n1, x_n1, y_n1 = t_n - h, x_n - h, y_n - h
    t_n2, x_n2, y_n2 = t_n - 2 * h, x_n - 2 * h, y_n - 2 * h
    t_n3, x_n3, y_n3 = t_n - 3 * h, x_n - 3 * h, y_n - 3 * h

    f_n = f_x(t_n, x_n, y_n)
    f_n1 = f_x(t_n1, x_n1, y_n1)
    f_n2 = f_x(t_n2, x_n2, y_n2)
    f_n3 = f_x(t_n3, x_n3, y_n3)

    g_n = f_y(t_n, x_n, y_n)
    g_n1 = f_y(t_n1, x_n1, y_n1)
    g_n2 = f_y(t_n2, x_n2, y_n2)
    g_n3 = f_y(t_n3, x_n3, y_n3)

    x_next = x_n + (h / 24) * (9 * f_n + 19 * f_n1 - 5 * f_n2 + f_n3)
    y_next = y_n + (h / 24) * (9 * g_n + 19 * g_n1 - 5 * g_n2 + g_n3)

    return x_next, y_next

# Метод Runge-Kutta 3-го порядку
def runge_kutta_3(t_n, x_n, y_n, h):
    K1_x = f_x(t_n, x_n, y_n)
    K1_y = f_y(t_n, x_n, y_n)

    K2_x = f_x(t_n + h / 3, x_n + h / 3 * K1_x, y_n + h / 3 * K1_y)
    K2_y = f_y(t_n + h / 3, x_n + h / 3 * K1_x, y_n + h / 3 * K1_y)

    K3_x = f_x(t_n + 2 * h / 3, x_n + 2 * h / 3 * K2_x, y_n + 2 * h / 3 * K2_y)
    K3_y = f_y(t_n + 2 * h / 3, x_n + 2 * h / 3 * K2_x, y_n + 2 * h / 3 * K2_y)

    x_next = x_n + (h / 4) * (K1_x + 3 * K3_x)
    y_next = y_n + (h / 4) * (K1_y + 3 * K3_y)

    return x_next, y_next

# Основна функція для розв'язку
def solve_differential_equations():
    k = 1000
    t0 = 0
    t1 = 5 / 1998  # Оцінка часу t1
    
    x_0 = 1
    y_0 = 1
    
    h = t1/k  # Початковий крок
    epsilon = 0.000001  # Точність
    
    # Початкові значення
    t_values = [t0]
    x_values = [x_0]
    y_values = [y_0]

    # Перехід за допомогою Runge-Kutta
    for _ in range(3):
        x_next, y_next = runge_kutta_3(t_values[-1], x_values[-1], y_values[-1], h)
        t_values.append(t_values[-1] + h)
        x_values.append(x_next)
        y_values.append(y_next)

    # Основний цикл Adams-Bashforth 4-го порядку
    while t_values[-1] < t1:
        x_next, y_next = adams_bashforth_4(t_values[-1], x_values[-1], y_values[-1], h)
        t_values.append(t_values[-1] + h)
        x_values.append(x_next)
        y_values.append(y_next)

        # Перевірка збіжності
        if (abs(x_next - x_values[-2]) + abs(y_next - y_values[-2])) < epsilon:
            break
    
    # Зміна кроку h для наступного інтервалу
    h = (t1 - t0) / 3
    T = 1  # Загальний інтервал

    # Подальший перехід з новим h
    while t_values[-1] < T:
        x_next, y_next = adams_bashforth_4(t_values[-1], x_values[-1], y_values[-1], h)
        t_values.append(t_values[-1] + h)
        x_values.append(x_next)
        y_values.append(y_next)

        # Перевірка збіжності
        if (abs(x_next - x_values[-2]) + abs(y_next - y_values[-2])) < epsilon:
            break

    return t_values, x_values, y_values

# Запуск розв'язку
t_values, x_values, y_values = solve_differential_equations()

# Друк результатів
print("RESULT:")
for t, x, y in zip(t_values, x_values, y_values):
    exact_x = 4*math.exp(-t) - 3*math.exp(-1000*t)
    exact_y = 2*math.exp(-t) + 3*math.exp(-1000*t)
    print(f"t = {t:.6f}, x = {x:.6f}, y = {y:.6f}")
    print(f"Точність: |x(t) - exact_x| = {abs(x - exact_x):.6f}, |y(t) - exact_y| = {abs(y - exact_y):.6f}")
