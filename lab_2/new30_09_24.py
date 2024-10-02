import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Визначення функцій
def f_1(Q):
    return 24 * (np.cos(2 * Q) - np.cos(3 * Q))

def f_2(Q):
    return 24 * (np.sin(2 * Q) - np.sin(3 * Q))

def f_3(Q):
    return 9 * np.cos(3 * Q) + 19 * np.cos(2 * Q) - 5 * np.cos(Q) + 1

def f_4(Q):
    return 9 * np.sin(3 * Q) + 19 * np.sin(2 * Q) - 5 * np.sin(Q)

def U(Q):
    return (f_1(Q) * f_3(Q) + f_2(Q) * f_4(Q)) / (f_3(Q)**2 + f_4(Q)**2)

def V(Q):
    return (f_2(Q) * f_3(Q) + f_1(Q) * f_4(Q)) / (f_3(Q)**2 + f_4(Q)**2)

# Параметри
h = 0.1  # Крок
N = 2000   # Кількість кроків
y = np.zeros(N)  # Масив для зберігання значень y
Q_values = np.linspace(0, 2 * np.pi, N)

# Розрахунок y[n] за формулою НА4
for n in range(3, N):
    f_n = f_3(Q_values[n])
    f_n1 = f_3(Q_values[n-1])
    f_n2 = f_3(Q_values[n-2])
    f_n3 = f_3(Q_values[n-3])
    
    y[n] = y[n-1] + (h/24) * (9*f_n3 + 19*f_n2 - 5*f_n1 + f_n)

# Побудова таблиці значень
U_values = U(Q_values)
V_values = V(Q_values)
data = pd.DataFrame({'Q': Q_values, 'U': U_values, 'V': V_values, 'y': y})
print(data)

# Побудова графіку V = f(U)
plt.figure(figsize=(10, 6))
plt.plot(U_values, V_values, label='V = f(U)', color='blue', marker='o')
plt.axhline(1, color='red', linestyle='--', label='V = 1')
plt.axhline(-1, color='red', linestyle='--', label='V = -1')
plt.title('Графік V = f(U)')
plt.xlabel('U')
plt.ylabel('V')
plt.legend()
plt.grid()
plt.show()

# Визначення області абсолютної стійкості
stable_region = data[(data['V'] > -1) & (data['V'] < 1)]
print("Область абсолютної стійкості (Q, U, V, y):")
print(stable_region)
