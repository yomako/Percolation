import codecs
import time
import random
from calculations import average_degree, clusters_identifier
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def complexity_power_function(x, C, alpha):
    """
    One-dimensional power function.

    Args:
        alpha: index
        C: normalization factor
        x: function's argument
    """
    return C * np.power(x, alpha)


def calculate_complexity():
    """Calculates computational complexity and puts results to text file."""
    size_array = []
    elapsed_time = []
    file = codecs.open("data/complexity.txt", "w", "utf-8")

    for s in range(10, 500):
        size = s*10
        size_array.append(size)
        hold = 0.001
        neighbourhood_matrix = np.zeros(shape=(size, size))
        for i in range(size):
            for j in range(size - 1 - i):
                rand = random.uniform(0, 1)
                if rand < hold:
                    neighbourhood_matrix[i][size - 1 - j] = 1
                    neighbourhood_matrix[size - 1 - j][i] = 1
                else:
                    neighbourhood_matrix[i][size - 1 - j] = 0
                    neighbourhood_matrix[size - 1 - j][i] = 0

        x = []
        for i in range(size):
            x.append(i)

        t = time.time()
        clusters_identifier(neighbourhood_matrix)
        print(size, average_degree(neighbourhood_matrix))
        elapsed = time.time() - t
        elapsed_time.append(elapsed)
        file.write("{0};{1};\n".format(size, elapsed))


def show_computational_complexity_diagram():
    """Shows computational complexity diagram."""
    x = []
    y = []
    for line3 in enumerate(codecs.open('data/complexity.txt', "r", "utf-8")):
        splited = line3[1].split(';')
        y.append(float(splited[1]))
        x.append(float(splited[0]))
    popt, pcov = curve_fit(complexity_power_function, x, y)
    print("C", popt[0], np.sqrt(pcov[0][0]))
    print("alfa", popt[1], np.sqrt(pcov[1][1]))
    chi_y = complexity_power_function(x, *popt)
    chi_square = 0
    for i in range(len(y)):
        chi_square += (y[i] - chi_y[i]) ** 2 / chi_y[i]
    print("chi^2", chi_square)
    print("doF", len(y)-1)
    print("chi^2/doF", chi_square/(len(y)-1))
    fit_x = np.linspace(100, x[len(x) - 1])
    plt.plot(x, y, '*', fit_x, complexity_power_function(fit_x, *popt), 'r-', label='fit')
    plt.title("Computational complexity as a function of size of the net")
    plt.xlabel("size N")
    plt.ylabel("Operating time of algorithm (s)")
    plt.show()
