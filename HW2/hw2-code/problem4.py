'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Probabilistic Graphical Models

Homework 2, Problem 4 - Inference in Ising Model

Haekyu Park
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''

import numpy as np


def main(N=10):

    # Initialize p
    p = np.zeros(2**N)
    for i in range(2**N):
        X = decimal_to_binary(i)
        p[i] = np.exp(np.sum(X[:-1] * X[1:]))

    # Initialize phi
    phi = np.zeros((2**N, 2**N))
    for i in range(0, 2**N):
        X_next = decimal_to_binary(i)
        for j in range(0, 2**N):
            X = decimal_to_binary(j)
            if len(X_next) != len(X):
                phi[i, j] = np.exp(0)
            else:
                phi[i, j] = np.exp(np.sum(X == X_next))

    # Message passing
    prev_msg = np.ones(2**N)
    curr_msg = np.zeros(2**N)

    for i in range(N - 1):

        curr_msg = np.zeros(2**N)

        for j in range(2**N):
            val = 0
            for k in range(2**N):
                val += p[k] * phi[j, k] * prev_msg[k]
            curr_msg[j] = val

        # Update prev_msg
        prev_msg = np.copy(curr_msg)

    # Get the normalization Z
    Z = 0
    for j in range(2**N):
        Z += p[j] * curr_msg[j]
    Z = np.log(Z)

    print(Z)


def decimal_to_binary(n):

    if n == 0:
        return np.array([0])

    elif n == 1:
        return np.array([1])

    else:
        result_before = decimal_to_binary(n // 2)
        return np.concatenate((result_before, np.array([n % 2])))


if __name__ == "__main__":
    N = 10
    main(N)
