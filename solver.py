import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


def upper_diagonal(h, i):
    result = -1 / h
    x_k = i * h
    print("Licze UD dla x_k: {}".format(x_k))
    if x_k + h <= 1:
        print("Całość k(x) = 1")
        result += 0.5
    elif x_k >= 1:
        print("Całość k(x) = 2")
        result += 1
    else:
        print("Punkt a = 1 między {} i {}".format(x_k, x_k + h))
        result += 5 / 8
    print("Rezultat: {}".format(result))
    return result


def lower_diagonal(h, i):
    result = -1 / h
    x_k = i * h
    print("Licze LD dla x_k: {}".format(x_k))
    if x_k <= 1:
        print("Całość k(x) = 1")
        result -= 0.5
    elif x_k - h >= 1:
        print("Całość k(x) = 2")
        result -= 1
    else:
        print("Punkt a = 1 między {} i {}".format(x_k - h, x_k))
        result -= 7 / 8
    print("Rezultat: {}".format(result))
    return result


def main_diagonal(h, i):
    result = 2 / h
    x_k = i * h
    print("Licze MD dla x_k: {}".format(x_k))
    if x_k == 1:
        print("x_k == 1")
        result -= 0.5
    elif ((x_k - h < 1) and (x_k > 1)) or ((x_k < 1) and (x_k + h > 1)):
        print("Przypadek x_k-1: {} x_k: {} x_k+1: {} a = 1".format(x_k - h, x_k, x_k + h))
        result -= 1 / 8
    print("Rezultat: {}".format(result))
    return result


def mes(g0, u2, f, n):
    h = 2.0 / n
    print("h: {}".format(h))
    matrix = np.zeros((n, n))
    for i in range(n):
        matrix[i][i] = main_diagonal(h, i)
        if i > 0:
            matrix[i - 1][i] = upper_diagonal(h, i - 1)
            matrix[i][i - 1] = lower_diagonal(h, i)
    matrix[0][0] = 1 / h - 0.5
    vector_r = np.zeros((n, 1))
    for i in range(n):
        vector_r[i] = 4 / 3 * f(h * i) * h
    vector_r[0] = (f(0) + 2 * f(h / 2)) * h / 6
    vector_r[0] -= g0
    vector_r[n - 1] -= u2 * (-1 / h + 1)
    result = la.solve(matrix, vector_r)
    print(matrix)
    print(vector_r)
    print(result)
    points = np.linspace(0.0, 2.0, n + 1)
    values = np.zeros(n + 1)
    for i in range(n):
        values[i] = result[i]
    values[n] = u2
    plt.plot(points, values)
    plt.show()


def main():
    mes(0, -3, lambda x: x ** 2 + 10 * x + np.exp(x), 1002)


if __name__ == "__main__":
    main()
