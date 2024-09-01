import numpy as np

# Задаємо функції x' і y'
def f_x(t, x, y):
    return t / y

def f_y(t, x, y):
    return -t / x

# Метод Ньютона для уточнення розв'язку неявної схеми
def newton_method(F, dF, initial_guess, tol=1e-10, max_iter=100):
    x = initial_guess
    for _ in range(max_iter):
        fx = F(x)
        dfx = dF(x)
        dx = np.linalg.solve(dfx, -fx)
        x = x + dx
        if np.linalg.norm(dx) < tol:
            break
    return x

# Функція для обчислення за схемою ЯРК3
def RK3_step(t, x, y, h):
    K1_x = f_x(t, x, y)
    K1_y = f_y(t, x, y)
    
    K2_x = f_x(t + h / 3, x + (h / 3) * K1_x, y + (h / 3) * K1_y)
    K2_y = f_y(t + h / 3, x + (h / 3) * K1_x, y + (h / 3) * K1_y)
    
    K3_x = f_x(t + 2 * h / 3, x + (2 / 3) * h * K2_x, y + (2 / 3) * h * K2_y)
    K3_y = f_y(t + 2 * h / 3, x + (2 / 3) * h * K2_x, y + (2 / 3) * h * K2_y)
    
    x_next = x + (h / 4) * (K1_x + 3 * K3_x)
    y_next = y + (h / 4) * (K1_y + 3 * K3_y)
    
    return x_next, y_next

# Основна програма для розв'язання задачі
def solve_na4_rk3(t0, x0, y0, h, N):
    t_values = [t0]
    x_values = [x0]
    y_values = [y0]
    
    # Ініціалізуємо початкові значення за допомогою схеми ЯРК3
    for _ in range(3):
        x_next, y_next = RK3_step(t_values[-1], x_values[-1], y_values[-1], h)
        t_values.append(t_values[-1] + h)
        x_values.append(x_next)
        y_values.append(y_next)
    
    # Головна схема НА4
    for n in range(3, N):
        t_n3 = t_values[n - 3]
        t_n2 = t_values[n - 2]
        t_n1 = t_values[n - 1]
        t_n = t_values[n]
        
        x_n3, x_n2, x_n1, x_n = x_values[n - 3], x_values[n - 2], x_values[n - 1], x_values[n]
        y_n3, y_n2, y_n1, y_n = y_values[n - 3], y_values[n - 2], y_values[n - 1], y_values[n]
        
        def F(z):
            x_next, y_next = z
            fx_next = f_x(t_n + h, x_next, y_next)
            fy_next = f_y(t_n + h, x_next, y_next)
            
            return np.array([
                x_next - x_n2 - (h / 24) * (9 * fx_next + 19 * f_x(t_n, x_n, y_n) - 5 * f_x(t_n1, x_n1, y_n1) + f_x(t_n2, x_n2, y_n2)),
                y_next - y_n2 - (h / 24) * (9 * fy_next + 19 * f_y(t_n, x_n, y_n) - 5 * f_y(t_n1, x_n1, y_n1) + f_y(t_n2, x_n2, y_n2))
            ])
        
        def dF(z):
            x_next, y_next = z
            fx_next = f_x(t_n + h, x_next, y_next)
            fy_next = f_y(t_n + h, x_next, y_next)
            
            dfx_dx = 0
            dfx_dy = -t_n / y_next**2
            dfy_dx = t_n / x_next**2
            dfy_dy = 0
            
            return np.array([
                [1 - (h / 24) * 9 * dfx_dx, -(h / 24) * 9 * dfx_dy],
                [-(h / 24) * 9 * dfy_dx, 1 - (h / 24) * 9 * dfy_dy]
            ])
        
        # Використовуємо метод Ньютона для уточнення розв'язку
        initial_guess = np.array([x_n, y_n])
        x_next, y_next = newton_method(F, dF, initial_guess)
        
        t_values.append(t_n + h)
        x_values.append(x_next)
        y_values.append(y_next)
    
    return t_values, x_values, y_values

# Початкові умови
t0 = 0
x0 = 1
y0 = 1
h = 0.1
N = 1000


# Розв'язуємо задачу
t_values, x_values, y_values = solve_na4_rk3(t0, x0, y0, h, N)

# Виводимо результати
for t, x, y in zip(t_values, x_values, y_values):
    x_exact = np.exp(t**2 / 2)
    y_exact = np.exp(-t**2 / 2)
    print(f"t = {t:.2f}, x = {x:.5f}, x_exact = {x_exact:.5f}, y = {y:.5f}, y_exact = {y_exact:.5f}")

    if t > 1:
        break
