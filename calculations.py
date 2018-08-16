import numpy as np
import random
import codecs
import math
import time


class Cluster:
    """
    Symbolizes cluster in the net.

    Attributes:
        nodes: list of nodes within the cluster
    """
    nodes = None

    def __init__(self):
        """Inits Cluster"""
        self.nodes = []

    def add_to_cluster(self, node):
        """
        Adds node to nodes list.

        Args:
            node: index of node
        """
        self.nodes.append(node)

    def length(self):
        """Returns length of nodes list."""
        return len(self.nodes)

    def is_in(self, j):
        """
        Checks if node with id j is in nodes list.

        Args:
            j: index of cluster

        Returns:
            True if node with id j is in nodes list, else returns False.
        """
        if j in self.nodes:
            return True
        else:
            return False

    def get_nodes(self):
        """Returns nodes list."""
        return self.nodes

    def print_nodes(self):
        """Prints nodes list."""
        print(self.nodes)


def test(clusters, size):
    """
    Checks if sum of clusters nodes list is equal to size parameter.

    Args:
        clusters: list of clusters
        size: size of the net

    Returns:
        True if sum of clusters nodes list is equal to size parameter, else returns False.
    """
    suma = 0
    for b in clusters:
        suma += b.length()
    if suma == size:
        return True
    else:
        return False


def nodes_degrees(matrix):
    """
    Calculates sum of degrees in the matrix.

    Args:
        matrix: matrix of neighborhood of analysed net

    Returns:
        Sum of degrees in the matrix.
    """
    return [int(sum(matrix[i])) for i in range(matrix.shape[0])]


def average_degree(matrix):
    """
    Calculates average degree ot analysed net.

    Args:
        matrix: matrix of neighborhood of analysed net

    Returns:
         Average degree of the net."""
    degrees = [int(sum(matrix[i])) for i in range(matrix.shape[0])]
    return sum(degrees)/matrix.shape[0]


def the_biggest_cluster_size(clusters):
    """
    Calculates size of the biggest cluster n analysed net.

    Args:
        clusters: list of clusters

    Returns:
         Size of the biggest cluster."""
    temp = []
    for b in clusters:
        temp.append(b.length())
    temp.sort(reverse=True)
    return temp[0]


def average_cluster_size(clusters):
    """
    Calculates average cluster size.

    Args:
        clusters: list of clusters

    Returns:
         Average cluster size.
    """
    temp = []
    for cl in clusters:
        temp.append(cl.length())
    temp.sort(reverse=True)
    temp.pop(0)
    if len(temp) == 0:
        return 0
    else:
        return sum(temp)/len(temp)


def normalize(array, size):
    """
    Normalizes array.

    Args:
        array: array to normalization
        size: size of neighborhood matrix

    Returns:
         Normalized array.
    """
    normalized_array = []
    for s in array:
        normalized_array.append(s/size)
    return normalized_array


def make_neighborhood_matrix(size, hold):
    """
    Makes random neighborhood matrix.

    Args:
        size: size of neighborhood matrix
        hold: probability of drawing edge between two nodes

    Returns:
        Neighborhood matrix.
    """
    neighborhood_matrix = np.zeros(shape=(size, size))
    for i in range(size):
        for j in range(size - 1 - i):
            rand = random.uniform(0, 1)
            if rand < hold:
                neighborhood_matrix[i][size - 1 - j] = 1
                neighborhood_matrix[size - 1 - j][i] = 1
            else:
                neighborhood_matrix[i][size - 1 - j] = 0
                neighborhood_matrix[size - 1 - j][i] = 0
    return neighborhood_matrix


def clusters_identifier(neighborhood_matrix):
    """
    Algorithm for identifying clusters in analysed net.

    Args:
         neighborhood_matrix: matrix of neighborhood of analysed net

    Returns:
        List of clusters.
    """
    clusters = []
    x = []
    for i in range(neighborhood_matrix.shape[0]):
        x.append(i)
    while True:
        try:
            pointer = x[0]
            temp_cluster = Cluster()
            temp_cluster.add_to_cluster(pointer)
            x.remove(pointer)

            cluster_neighbors = []

            for j in range(len(neighborhood_matrix[pointer])):
                if neighborhood_matrix[pointer][j] == 1 and temp_cluster.is_in(j) is False:
                    cluster_neighbors.append(j)
                    try:
                        x.remove(j)
                    except Exception:
                        pass

            h = len(cluster_neighbors)

            while h > 0:
                for i in cluster_neighbors:
                    temp_cluster.add_to_cluster(i)
                cluster_neighbors.clear()

                for i in range(h):
                    for j in range(len(neighborhood_matrix[int(temp_cluster.nodes[len(temp_cluster.nodes) - 1 - i])])):
                        if int(neighborhood_matrix[int(temp_cluster.nodes[len(temp_cluster.nodes) - 1 - i])][
                                   j]) == 1 and temp_cluster.is_in(j) is False and (j in cluster_neighbors) is False:
                            cluster_neighbors.append(j)
                            try:
                                x.remove(j)
                            except Exception:
                                pass

                h = len(cluster_neighbors)

            clusters.append(temp_cluster)
        except Exception:
            break
    return clusters


def calculations(events, size, hold=0.00001):
    """
    Main function which analyzes neighborhood matrix and puts results to text files.

    Args:
        events: number of calculative events
        size: size of the net
        hold: probability of drawing edge between two nodes

    """
    for q in range(events):
        step = math.pow(10, (-1*math.log10(size)-1.69897))
        file = codecs.open("data/data_{0}.txt".format(size), "a", "utf-8")
        file2 = codecs.open("data/clusters_histogram_{0}.txt".format(size), "a", "utf-8")
        file3 = codecs.open("data/average_clusters_{0}.txt".format(size), "a", "utf-8")
        the_biggest_clusters_sizes = []
        average_degrees = []
        elapsed_time = []
        for h in range(200):
            print(h)
            neighborhood_matrix = make_neighborhood_matrix(size, hold)
            t = time.time()
            clusters = clusters_identifier(neighborhood_matrix)
            elapsed = time.time() - t
            elapsed_time.append(elapsed)
            the_biggest_clusters_sizes.append(the_biggest_cluster_size(clusters))
            average_degrees.append(average_degree(neighborhood_matrix))
            line = "{0};{1};\n".format(average_degree(neighborhood_matrix), the_biggest_cluster_size(clusters)/size)
            line3 = "{0};{1};\n".format(average_degree(neighborhood_matrix), average_cluster_size(clusters))
            if 0.99 <= average_degree(neighborhood_matrix) <= 1.01:
                for cluster in clusters:
                    line2 = "{0};\n".format(cluster.length())
                    file2.write(line2)

            file.write(line)
            file3.write(line3)
            hold += step

        the_biggest_clusters_sizes = normalize(the_biggest_clusters_sizes, size)
        print(len(the_biggest_clusters_sizes))
        the_biggest_clusters_sizes.sort()
        average_degrees.sort()
        print("elapsed", sum(elapsed_time)/len(elapsed_time))
