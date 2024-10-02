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


# Побудова таблиці значень
U_values = U(Q_values)
V_values = V(Q_values)
data = pd.DataFrame({'Q': Q_values, 'U': U_values, 'V': V_values, 'y': y})
print(data)
p = 1/2
sigma = (24*(p*p - p*p*p))/(9*p*p*p + 19*p*p - 5*p + 1)

# Побудова графіку V = f(U)
plt.figure(figsize=(10, 6))
plt.plot(U_values, V_values, label='V = f(U)', color='blue', marker='o')
plt.plot(sigma, 0, marker='o', label='sigma', color='r')  
plt.fill_betweenx(np.linspace(-5, 5, 400), -5, 5, color='blue', alpha=0.2)  # Заповнюємо фон
plt.fill(U_values, V_values, color='white')
plt.title('Графік V = f(U)')
plt.xlabel('U')
plt.ylabel('V')
plt.legend()
plt.grid()
plt.show()



