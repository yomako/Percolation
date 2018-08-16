import codecs
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit


def power_function(x, C, alpha):
    """
    One-dimensional power function.

    Args:
        alpha: index
        C: normalization factor
    """
    return C * np.power(x, -alpha)


def load_percolation_data(size, bins=200, tbcs_n_returns=True):
    """
    Loads data of percolational phase transition from data_{size}.txt file.

    Args:
        size: size of considered net
        bins: number of histogram's bins
        tbcs_n_returns: if True returns table of the biggest clusters sizes

    Returns:
        If tbcs_n_returns is True returns table of degrees, corresponding the biggest clusters sizes values and errors,
            if not returns table of degrees and table of the biggest clusters sizes errors only.
    """
    the_biggest_clusters_sizes = []
    average_degrees = []
    degrees = []
    for line in enumerate(codecs.open('data/data_{0}.txt'.format(size), "r", "utf-8")):
        divided = line[1].split(';')
        the_biggest_clusters_sizes.append(float(divided[1]))
        average_degrees.append(float(divided[0]))
    tbcs_n = []
    tbcs_err = []
    low = min(average_degrees)
    for i in range(bins):
        high = low + (max(average_degrees) - min(average_degrees)) / 200
        bin = []
        for j in range(len(average_degrees)):
            if low <= average_degrees[j] < high:
                bin.append(the_biggest_clusters_sizes[j])
        suma = 0
        if len(bin) > 1:
            mean = sum(bin) / len(bin)
            for j in bin:
                suma += abs(mean - j) * abs(mean - j)
            suma /= (len(bin) - 1)
            suma = math.sqrt(suma)
            tbcs_n.append(mean)
            tbcs_err.append(suma)
            degrees.append((high + low) / 2)

        low = high

    if tbcs_n_returns is True:
        return degrees, tbcs_n, tbcs_err
    else:
        return degrees, tbcs_err


def show_percolation_diagram(size, bins, start):
    """
    Shows diagram of percolational phase transition - size of the biggest cluster as a function of average degree.

    Args:
        size: size of considered net
        bins: number of histogram's bins
        start: optional variable to figures numbering
    """
    degrees, tbcs_n, tbcs_err = load_percolation_data(size, bins)
    plt.figure(num=start + 1, figsize=(8, 6), dpi=80)
    plt.errorbar(degrees, tbcs_n, fmt='r', label="data", xerr=0, yerr=tbcs_err, ecolor='black', )
    plt.title("Diagram of percolational phase transition for net with number of nodes N={0}".format(size))
    plt.xlabel("<k>")
    plt.ylabel("$n_{S}$")
    plt.axis((min(degrees), max(degrees), min(tbcs_n), max(tbcs_n) * 1.1))


def show_diagram_of_fluctuations(size, bins, start):
    """
    Shows diagram of fluctuations of size of the biggest cluster as a function of average degree.

    Args:
        size: size of considered net
        bins: number of histogram's bins
        start: optional variable to figures numbering
    """
    degrees, tbcs_err = load_percolation_data(size, bins, False)
    plt.figure(num=start + 2, figsize=(8, 6), dpi=80)
    plt.plot(degrees, tbcs_err)
    plt.title("Diagram of fluctuations for net with number of nodes N={0}".format(size))
    plt.xlabel("<k>")
    plt.ylabel("$\Delta$S")
    plt.axis((min(degrees), max(degrees), min(tbcs_err), max(tbcs_err) * 1.1))


def load_clusters_size_data(size):
    """
    Loads data of clusters sizes distribution in phase transition point from clusters_histogram_{size}.txt file.

    Args:
        size: size of considered net

    Returns:
        List of clusters sizes and list of corresponding probabilities of finding that clusters in the net.
    """
    clusters_sizes_raw = []
    for line2 in enumerate(codecs.open('data/clusters_histogram_{0}.txt'.format(size), "r", "utf-8")):
        divided = line2[1].split(';')
        clusters_sizes_raw.append(int(float(divided[0])))

    clusters_sizes = []
    probability_of_cluster_sizes = []
    for i in range(1, 8):
        ii = 0
        clusters_sizes.append(i)
        for k in clusters_sizes_raw:
            if k == i:
                ii += 1
        probability_of_cluster_sizes.append(ii)

    norm = sum(clusters_sizes_raw)
    for i in range(len(probability_of_cluster_sizes)):
        probability_of_cluster_sizes[i] /= norm

    return clusters_sizes, probability_of_cluster_sizes


def fit_to_clusters_size_distribution(clusters_sizes, probability_of_cluster_sizes, print_parameters):
    """
    Fits power function to clusters size distribution.

    Args:
        clusters_sizes: list of clusters sizes
        probability_of_cluster_sizes: list of corresponding probabilities of finding that clusters in the net
        print_parameters: if is True prints chi^2 and values with errors of fitted C and alpha parameters

    Returns:
        Array of fitted values and list of fitted parameters.
    """
    print(clusters_sizes)
    print(probability_of_cluster_sizes)
    popt, pcov = curve_fit(power_function, clusters_sizes, probability_of_cluster_sizes)
    fit_x = np.linspace(1, clusters_sizes[len(clusters_sizes) - 1])
    chi_y = power_function(clusters_sizes, *popt)
    chi_square = 0
    for i in range(len(clusters_sizes)):
        chi_square += (probability_of_cluster_sizes[i] - chi_y[i]) ** 2 / chi_y[i]

    if print_parameters is True:
        print("T", chi_square)
        print("C = {0}({1})".format(popt[0], np.sqrt(pcov[0][0])))
        print("alpha = {0}({1})".format(popt[1], np.sqrt(pcov[1][1])))

    return fit_x, popt


def show_clusters_size_distribution(size, print_parameters, start):
    """
    Shows diagram of clusters size distribution.

    Args:
        size: size of considered net
        print_parameters: if is True prints chi^2 and values with errors of fitted C and alpha parameters
        start: optional variable to figures numbering
    """
    plt.figure(num=start + 3, figsize=(8, 6), dpi=80)
    clusters_sizes, probability_of_cluster_sizes = load_clusters_size_data(size)
    fit_x, popt = fit_to_clusters_size_distribution(clusters_sizes, probability_of_cluster_sizes, print_parameters)
    plt.plot(clusters_sizes, probability_of_cluster_sizes, '*', fit_x, power_function(fit_x, *popt), 'r-')
    plt.title("Average cluster size distribution for net with N={0}".format(size))
    plt.xlabel("S")
    plt.ylabel("p(S)")
    plt.legend(('Results of numerical calculations', 'Fitted power distribution'))
    plt.axis((min(clusters_sizes) - 0.5, max(clusters_sizes), min(probability_of_cluster_sizes),
              max(probability_of_cluster_sizes) + 0.05))


def load_average_cluster_size_data(size, bins=40):
    """
    Loads data of average clusters sizes as a function of average degree from average_clusters_{size}.txt file.

    Args:
        size: size of considered net
        bins: number of histogram's bins

    Returns:
        List of average degrees and list of corresponding average clusters sizes.
    """
    average_degrees_raw = []
    average_clusters_sizes_raw = []
    for line3 in enumerate(codecs.open('data/average_clusters_{0}.txt'.format(size), "r", "utf-8")):
        divided = line3[1].split(';')
        average_clusters_sizes_raw.append(float(divided[1]))
        average_degrees_raw.append(float(divided[0]))

    average_degrees = []
    average_clusters_sizes = []

    low = min(average_degrees_raw)
    for i in range(bins):
        high = low + (max(average_degrees_raw) - min(average_degrees_raw)) / bins
        bin = []
        for j in range(len(average_degrees_raw)):
            if low <= average_degrees_raw[j] < high and average_clusters_sizes_raw[j] >= 1:
                bin.append(average_clusters_sizes_raw[j])
        if bin:
            mean = sum(bin) / len(bin)
            average_clusters_sizes.append(mean)
            average_degrees.append((high + low) / 2)

        low = high

    return average_degrees, average_clusters_sizes


def show_average_cluster_size_data(size, start):
    """
    Shows diagram of average clusters sizes as a function of average degree.

    Args:
        size: size of considered net
        start: optional variable to figures numbering
    """
    plt.figure(num=start + 4, figsize=(8, 6), dpi=80)
    average_degree, average_cluster_size = load_average_cluster_size_data(size)
    plt.plot(average_degree, average_cluster_size)
    plt.title("Average cluster size for net with N={0}".format(size))
    plt.xlabel("<k>")
    plt.ylabel("<S>")
    plt.axis((min(average_degree), max(average_degree), 1, max(average_cluster_size) + 0.1))


def show_diagrams(size, print_parameters=True, start=0, percolation_diagram_bins=200, fluctuations_diagram_bins=200):
    """
    Shows all diagrams.

    Args:
        size: size of considered net
        print_parameters: if is True prints chi^2 and values with errors of fitted C and alpha parameters
        start: optional variable to figures numbering
        percolation_diagram_bins: number of percolation diagram bins
        fluctuations_diagram_bins: number of fluctuations diagram bins
    """
    try:
        show_percolation_diagram(size, percolation_diagram_bins, start)
        show_diagram_of_fluctuations(size, fluctuations_diagram_bins, start)
        show_clusters_size_distribution(size, print_parameters, start)
        show_average_cluster_size_data(size, start)
        plt.show()
    except FileNotFoundError:
        print("File not found")
