import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('results.csv', delimiter=',', skip_header=1)

N = data[:, 0]
h = data[:, 1]
C_err = data[:, 2]
L2_err = data[:, 3]
prec_time = data[:, 4]
iter_time = data[:, 5]

# График сходимости
plt.figure(figsize=(8,6))
plt.loglog(h, C_err, 'o-', label='C-норма (max)')
plt.loglog(h, L2_err, 's-', label='L2-норма')
ref = L2_err[0] * (h / h[0])**2
plt.loglog(h, ref, 'k--', label='$O(h^2)$')
plt.xlabel('Шаг сетки $h$')
plt.ylabel('Ошибка')
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend()
plt.title('Сходимость метода конечных разностей')
plt.savefig('convergence_plot.png', bbox_inches='tight')
plt.show()

# График времени
plt.figure(figsize=(8,6))
dof = (N - 1)**2
plt.plot(dof, prec_time, 'o-', label='Время сборки ILU2')
plt.plot(dof, iter_time, 's-', label='Время итераций BiCGStab')
plt.xlabel('Размер СЛАУ')
plt.ylabel('Время (сек)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.title('Зависимость времени от размера задачи')
plt.savefig('time_plot.png', bbox_inches='tight')
plt.show()
