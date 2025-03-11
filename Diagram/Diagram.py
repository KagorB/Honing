
import matplotlib.pyplot as plt
# Чтение данных из файла
data_file = "../T_values.txt"
epsilon_star_values = []
T_values = []
A = 0  # Initialize A with a default value
with open(data_file, 'r') as file:
    lines = file.readlines()
    A = float(lines[0].strip())  # Read the value of A from the first line
    i = 1
    while i < len(lines):
        epsilon_star = float(lines[i].strip())
        epsilon_star_values.append(epsilon_star)
        T = list(map(float, lines[i + 1].split()))
        #T_values.append(sorted(T[-100:]))
        T_values.append(T[-100:])  # Используем последние 100 значений T
        i += 2

# Построение бифуркационной диаграммы
plt.figure(figsize=(10, 6))
for i in range(len(epsilon_star_values)):
    plt.scatter([epsilon_star_values[i]] * len(T_values[i]), T_values[i], s=1, color='blue')

plt.title(f"Бифуркационная диаграмма T от E* при Teta = {A}")
plt.xlabel("E*")
plt.ylabel("T")
plt.grid(True)
plt.show()