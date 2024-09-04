import math
import matplotlib.pyplot as plt

def f_x(t, x, y):
    return t / y

def f_y(t, x, y):
    return -t / x

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

def solve_differential_equations(t_start, t_end, x_0, y_0, h, epsilon):
    t_values = [t_start]
    x_values = [x_0]
    y_values = [y_0]
    
    x_1 = x_0 + h*f_x(t_start,x_0,y_0)
    y_1 = y_0 + h*f_y(t_start,x_0,y_0)
    
    while  True:
        x_11 = x_0 + (h/2)*f_x(t_start,x_0,y_0)
        y_11 = y_0 + (h/2)*f_y(t_start,x_0,y_0)
    
        x_12 = x_11 + (h/2)*f_x(t_start+(h/2),x_11,y_11)
        y_12 = y_11 + (h/2)*f_y(t_start+(h/2),x_11,y_11)
    
        if (abs(x_1 - x_12) + abs(y_1 -y_12) < epsilon):
            break
    
        h = (h/2)
        x_1 = x_11
        y_1 = y_11
        
    print(h)
    
     #Use Runge-Kutta 3 to calculate initial steps
    for _ in range(3):
        x_next, y_next = runge_kutta_3(t_values[-1], x_values[-1], y_values[-1], h)
        t_values.append(t_values[-1] + h)
        x_values.append(x_next)
        y_values.append(y_next)

    # Use Adams-Bashforth 4 to continue the solution
    while t_values[-1] < t_end:
        x_next, y_next = adams_bashforth_4(t_values[-1], x_values[-1], y_values[-1], h)
        t_values.append(t_values[-1] + h)
        x_values.append(x_next)
        y_values.append(y_next)

        # Check for convergence
        if abs(x_next - x_values[-2]) < epsilon and abs(y_next - y_values[-2]) < epsilon:
            break

    return t_values, x_values, y_values

# Test case
t_start = 0
t_end = 1
x_0 = 1
y_0 = 1
h = 0.3
epsilon = 0.01

t_values, x_values, y_values = solve_differential_equations(t_start, t_end, x_0, y_0, h, epsilon)
#Графік
x_points = []
y_points = []

x_points_test = []
y_points_test = []

# Друк результатів
print("")
print("RESULT:")
for t, x, y in zip(t_values, x_values, y_values):
    if t > 1:
        break
    exact_x = math.exp((t**2)/2)
    exact_y = math.exp(-((t**2)/2))
    
    x_points.append(x)
    y_points.append(y)
    
    x_points_test.append(exact_x)
    y_points_test.append(exact_y)
    
    print(f"t = {t:.4f}, x = {x:.4f}, y = {y:.4f},  exact_x = {abs(x - exact_x):.4f}, exact_y = {abs(y - exact_y):.4f}")

plt.plot(x_points, y_points, color="red")
plt.plot(x_points_test, y_points_test, color="blue")
plt.show()