import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin, pi

def G(p):
    return (24*(p*p*p - p*p))/(9*p*p*p + 19*p*p - 5*p + 1)

def f_1(Q):
    return 2*(cos(Q) - cos(2*Q))

def f_2(Q):
    return 2*(sin(Q) - sin(2*Q))

def f_3(Q):
    return 3*cos(Q) - 1

def f_4(Q):
    return sin(Q)

def U(Q):
    return (f_1(Q)*f_3(Q) + f_2(Q)*f_4(Q)) / (f_3(Q)**2 + f_4(Q)**2)

def V(Q):
    return (f_2(Q)*f_3(Q) + f_1(Q)*f_4(Q)) / (f_3(Q)**2 + f_4(Q)**2)

u_number = []
v_number = []
q_number = []
Q = 0
p = 4

while Q < 2*pi:
    q_number.append(Q)
    u_number.append(U(Q))
    v_number.append(V(Q))
    
    Q += 2*pi/100

p_points = G(p)
print(" ")
print(" p = ", p_points)
print(" ")

for q, u, v in zip(q_number, u_number, v_number):
    print(f" q = {q:.6f},  u = {u:.6f},  v = {v:.6f}")

# Створюємо графік
plt.figure(figsize=(8, 6))
plt.plot(u_number, v_number, label='Curve')
plt.scatter(G(p), 0, color='red', zorder=5, label=f'Point ({G(p)}, 0)')  # Точка
plt.xlabel('u')
plt.ylabel('v')
plt.title('Графік кривої з однією точкою')
plt.legend()
plt.grid(True)

# Показуємо графік
plt.show()
