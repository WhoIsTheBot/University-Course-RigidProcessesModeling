import numpy as np
import matplotlib.pyplot as plt

# Визначення функції для R(z)
def R(z, h):
    return 1 + (h / 24) * (9 * z**3 + 19 * z**2 - 5 * z + 1)

# Параметри
h = 1  # Крок
r_values = np.linspace(0, 2, 400)  # Дійсна частина
theta_values = np.linspace(0, 2 * np.pi, 400)  # Уявна частина
R_values = np.zeros((len(r_values), len(theta_values)), dtype=complex)

# Обчислення значень R(z) для комплексної площини
for i, r in enumerate(r_values):
    for j, theta in enumerate(theta_values):
        z = r * np.exp(1j * theta)
        R_values[i, j] = R(z, h)

# Створення маски для стійкої області
stability_region = np.abs(R_values) < 1

# Візуалізація
plt.figure(figsize=(8, 8))
plt.contourf(r_values[:, None] * np.cos(theta_values), r_values[:, None] * np.sin(theta_values), stability_region.T, levels=1, colors='lightblue')
plt.title('Область абсолютної стійкості НA 4')
plt.xlabel('Дійсна частина')
plt.ylabel('Уявна частина')
plt.axhline(0, color='black', lw=0.5)
plt.axvline(0, color='black', lw=0.5)
plt.grid()
plt.axis('equal')
plt.show()
