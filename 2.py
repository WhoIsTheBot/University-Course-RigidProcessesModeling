import math
import matplotlib.pyplot as plt

# Define the system of equations
def f_x(t, x, y):
    return 998 * x + 1998 * y

def f_y(t, x, y):
    return 999 * x - 1999 * y

# Adams-Bashforth Method (4th order)
def adams_bashforth_4(t_n, x_n, y_n, h):
    f_n = f_x(t_n, x_n, y_n)
    g_n = f_y(t_n, x_n, y_n)
    
    # Previous values (you would need to store these during iterations)
    f_n1 = f_x(t_n - h, x_n - h * f_n, y_n - h * g_n)
    g_n1 = f_y(t_n - h, x_n - h * f_n, y_n - h * g_n)

    # Continue for previous steps as needed...
    
    # Calculate next values
    x_next = x_n + (h / 24) * (9 * f_n + 19 * f_n1) # Simplified for demonstration
    y_next = y_n + (h / 24) * (9 * g_n + 19 * g_n1)

    return x_next, y_next

# Runge-Kutta Method (3rd order)
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

# Solving function
def solve_differential_equations(t_start, t_end, x_0, y_0, h):
    t_values = [t_start]
    x_values = [x_0]
    y_values = [y_0]

    # Initial Runge-Kutta steps
    for _ in range(3):
        x_next, y_next = runge_kutta_3(t_values[-1], x_values[-1], y_values[-1], h)
        t_values.append(t_values[-1] + h)
        x_values.append(x_next)
        y_values.append(y_next)

    # Use Adams-Bashforth to continue
    while t_values[-1] < t_end:
        x_next, y_next = adams_bashforth_4(t_values[-1], x_values[-1], y_values[-1], h)
        t_values.append(t_values[-1] + h)
        x_values.append(x_next)
        y_values.append(y_next)

    return t_values, x_values, y_values

# Test case parameters
t_start = 0
t_end = 5
x_0 = 1
y_0 = 1
h = 0.01 # Example step size

t_values, x_values, y_values = solve_differential_equations(t_start, t_end, x_0, y_0, h)

# Print results
print("\nRESULT:")
for t, x, y in zip(t_values, x_values, y_values):
    print(f"t = {t:.6f}, x = {x:.6f}, y = {y:.6f}")
