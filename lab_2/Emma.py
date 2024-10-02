import numpy as np
import matplotlib.pyplot as plt


 # Визначаємо параметри
t_values = np.linspace(0, 2 * np.pi, 200)  # 200 точок від 0 до 2π
u_values = ((24 * (np.cos(2 * t_values) - np.cos(3 * t_values))) * (9 * np.cos(3 * t_values) + 19 * np.cos(2 * t_values) - 5 * np.cos(t_values) + 1) + (24 * (np.sin(2 * t_values) - np.sin(3 * t_values))) * (9 * np.sin(3 * t_values) + 19 * np.sin(2 * t_values) - 5 * np.sin(t_values))) / ((9 * np.cos(3 * t_values) + 19 * np.cos(2 * t_values) - 5 * np.cos(t_values) + 1)**2 + (9 * np.sin(3 * t_values) + 19 * np.sin(2 * t_values) - 5 * np.sin(t_values))**2)
v_values = ((24 * (np.sin(2 * t_values) - np.sin(3 * t_values))) * (9 * np.cos(3 * t_values) + 19 * np.cos(2 * t_values) - 5 * np.cos(t_values) + 1) + (24 * (np.cos(2 * t_values) - np.cos(3 * t_values))) * (9 * np.sin(3 * t_values) + 19 * np.sin(2 * t_values) - 5 * np.sin(t_values))) / ((9 * np.cos(3 * t_values) + 19 * np.cos(2 * t_values) - 5 * np.cos(t_values) + 1)**2 + (9 * np.sin(3 * t_values) + 19 * np.sin(2 * t_values) - 5 * np.sin(t_values))**2)

p = 1/2
sigma = 1 - (18/(11*p) - 9/(11*p**2) + 2/(11*p**3))
print(sigma)

# Створюємо графік
plt.figure(figsize=(10, 6))

# Встановлюємо межі графіка
plt.xlim(-5, 5)
plt.ylim(-5, 5)

# Заповнюємо область поза фігурою
plt.fill_betweenx(np.linspace(-5, 5, 400), -5, 5, color='blue', alpha=0.2)  # Заповнюємо фон
plt.fill(u_values, v_values, color='white')  # Залишаємо фігуру білою

# З'єднуємо точки
plt.plot(u_values, v_values, marker='o', label='u vs v', color='b', markersize=1) 
plt.plot(sigma, 0, marker='o', label='sigma', color='r')   
plt.title('Графік залежності u(t) та v(t)')
plt.xlabel('u(t)')
plt.ylabel('v(t)')

# Додаємо осі
plt.axhline(0, color='black', lw=0.5, ls='--')  
plt.axvline(0, color='black', lw=0.5, ls='--')
plt.grid()
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')  # Встановлюємо рівні пропорції
plt.show()
