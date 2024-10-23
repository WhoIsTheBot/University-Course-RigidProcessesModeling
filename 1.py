import numpy as np

# Визначення системи
A = np.array([[998, 1998], [999, -1999]])
t0 = 0
t1 = 5
k = 50
h = t1 / k

# Початкові умови
x0 = 0  # Встановіть свої значення
y0 = 0

# Кількість кроків
steps = int((t1 - t0) / h)
x = np.zeros(steps)
y = np.zeros(steps)
x[0] = x0
y[0] = y0

# Головна схема НА4
for n in range(steps - 3):
    # Обчислюємо значення f
    f_n = A @ np.array([x[n], y[n]])
    f_n_plus_1 = A @ np.array([x[n + 1], y[n + 1]])
    f_n_plus_2 = A @ np.array([x[n + 2], y[n + 2]])
    f_n_plus_3 = A @ np.array([x[n + 3], y[n + 3]])

    # Обчислюємо нові значення x і y
    y[n + 3] = y[n + 2] + (h / 24) * (9 * f_n_plus_3[1] + 19 * f_n_plus_2[1] - 5 * f_n_plus_1[1] + f_n[1])
    x[n + 3] = x[n + 2] + (h / 24) * (9 * f_n_plus_3[0] + 19 * f_n_plus_2[0] - 5 * f_n_plus_1[0] + f_n[0])
    
    # Допоміжна ЯРК3
    K1 = f_n
    K2 = A @ (np.array([x[n] + (h / 3) * K1[0], y[n] + (h / 3) * K1[1]]))
    K3 = A @ (np.array([x[n] + (2/3) * h * K2[0], y[n] + (2/3) * h * K2[1]]))
    
    y[n + 1] = y[n] + (h / 4) * (K1[1] + 3 * K3[1])
    x[n + 1] = x[n] + (h / 4) * (K1[0] + 3 * K3[0])

# Виведення результатів
for i in range(steps):
    print(f't={t0 + i * h}, x={x[i]}, y={y[i]}')
