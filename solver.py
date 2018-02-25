import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


# Counting upper diagonal values
def upper_diagonal(h, i):
    result = -1 / h
    x_k = i * h
    if x_k + h <= 1:
        result += 0.5
    elif x_k >= 1:
        result += 1
    else:
        result += 5 / 8
    return result


# Counting lower diagonal values
def lower_diagonal(h, i):
    result = -1 / h
    x_k = i * h
    if x_k <= 1:
        result -= 0.5
    elif x_k - h >= 1:
        result -= 1
    else:
        result -= 7 / 8
    return result


# Counting main diagonal values
def main_diagonal(h, i):
    result = 2 / h
    x_k = i * h
    if x_k == 1:
        result -= 0.5
    elif ((x_k - h < 1) and (x_k > 1)) or ((x_k < 1) and (x_k + h > 1)):
        result -= 1 / 8
    return result


def fem(g0, u2, f, n):
    h = 2.0 / n
    matrix = np.zeros((n, n))
    for i in range(n):
        matrix[i][i] = main_diagonal(h, i)
        if i > 0:
            matrix[i - 1][i] = upper_diagonal(h, i - 1)
            matrix[i][i - 1] = lower_diagonal(h, i)
    matrix[0][0] = 1 / h - 0.5
    vector_r = np.zeros((n, 1))
    # Simpson approximation of integrals
    for i in range(n):
        vector_r[i] = 4 / 3 * f(h * i) * h
    # Special case B(e_0)
    vector_r[0] = (f(0) + 2 * f(h / 2)) * h / 6
    # Neumann condition
    vector_r[0] -= g0
    # Dirichlet shift
    vector_r[n - 1] -= u2 * (-1 / h + 1)
    result = la.solve(matrix, vector_r)

    # Showing results
    points = np.linspace(0.0, 2.0, n + 1)
    values = np.zeros(n + 1)
    for i in range(n):
        values[i] = result[i]
    values[n] = u2
    plt.plot(points, values)
    plt.show()


def main():
    fem(-10, -3, lambda x: np.log2(x + 1) + np.exp(x ** 3 - 17), 1112)


if __name__ == "__main__":
    main()
