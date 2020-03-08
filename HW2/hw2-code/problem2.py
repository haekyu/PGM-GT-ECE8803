'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Probabilistic Graphical Models

Homework 2 - Inference in HMM

Haekyu Park
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''


'''
Import packages
'''
import numpy as np
from scipy.io import loadmat


'''
Main
'''


def main():
    '''
    Load noisy string
    - V: visible sequence
    '''

    mat_file_path = './data/noisystring.mat'
    noisy_str = load_data(mat_file_path)
    V = visible_sequence(noisy_str)


    '''
    Initialize states
    - T: Transition matrix
    - E: Emission matrix
    '''

    # Firstnames and surnames
    firstnames = ['david', 'anton', 'fred', 'jim', 'barry']
    surnames = ['barber', 'ilsung', 'fox', 'chain', 'fitzwilliam', 'quinceadams', 'grafvonunterhosen']

    # We have 2 more states, for the start state and the middle state
    num_states = 2 + sum(map(lambda x: len(x), firstnames)) + sum(map(lambda x: len(x), surnames))

    # Initialize transition matrix T and emission matrix E
    num_alphabets = 26
    T, E = init_states(num_states, num_alphabets)

    # Firstname states
    state_th_2_firstname_th = {}
    state_th = firstname_states(firstnames, T, E, state_th_2_firstname_th)

    # Surname states
    state_th_2_surname_th = {}
    surname_states(surnames, T, E, state_th_2_surname_th, state_th)


    '''
    Initialize distributions
    - phi: distributions
    '''
    phi = np.zeros(num_states)
    phi[0] = 1

    '''
    Viterbi
    '''
    seq = viterbi(V, T, E, phi)

    '''
    Check occurrence
    '''
    M = build_occurrence_matrix(firstnames, surnames, seq, state_th_2_firstname_th, state_th_2_surname_th)
    print(M)



'''
Parse noisystring
'''



def load_data(mat_file_path):
    data = loadmat(mat_file_path)
    return data['noisystring'][0]


def char2num(c):
    return ord(c) - ord('a')


def visible_sequence(noisy_str):
    noisy_chars = list(noisy_str)
    return list(map(lambda x: char2num(x), noisy_chars))



'''
States
'''


def init_states(num_states, num_alphabets):

    T = np.zeros((num_states, num_states))
    E = np.ones((num_states, num_alphabets)) * 0.7 / (num_alphabets - 1)

    # The start state
    T[0, 0] = 0.8
    E[0, :] = 1 / num_alphabets

    # The middle state
    T[1, 1] = 0.8
    E[1, :] = 1 / num_alphabets

    return T, E


def firstname_states(firstnames, T, E, state_th_2_firstname_th):
    # Name state starts from 2 (0: start state, 1: middle state)
    state_th = 2

    for firstname_th, firstname in enumerate(firstnames):
        # We start to generate a firstname with probabilty 0.2
        T[0, state_th] = 0.2 * (1 / len(firstnames))
        state_th_2_firstname_th[state_th] = firstname_th

        # Update T, E
        for char_th, char in enumerate(firstname):

            # Update T
            if char_th < len(firstname) - 1:
                # During concatenating names
                T[state_th + char_th, state_th + char_th + 1] = 1
            else:
                # At the last char of the names
                T[state_th + char_th, 1] = 1

            # Update E
            E[state_th + char_th, char2num(char)] = 0.3

        # Update state_th
        state_th += len(firstname)

    return state_th


def surname_states(surnames, T, E, state_th_2_surname_th, state_th):
    for surname_th, surname in enumerate(surnames):
        # We start to generate a firstname with probabilty 0.2
        T[1, state_th] = 0.2 * 1 / len(surname)
        state_th_2_surname_th[state_th] = surname_th

        # Update T, E
        for char_th, char in enumerate(surname):

            # Update T
            if char_th < len(surname) - 1:
                # During concatenating names
                T[state_th + char_th, state_th + char_th + 1] = 1
            else:
                # At the last char of the names
                T[state_th + char_th, 0] = 1

            # Update E
            E[state_th + char_th, char2num(char)] = 0.3

        # Update state_th
        state_th += len(surname)


'''
Viterbi
'''


def viterbi(V, T, E, phi):

    # Get the shapes
    len_seqs = len(V)
    num_states = T.shape[0]

    # Initialize the likelihood
    omega = np.zeros((len_seqs, num_states))
    omega[0, :] = np.log(phi) + np.log(E[:, V[0]])

    # Initialize the probable state
    probable_states = np.ones((len_seqs, num_states))

    # Forward
    for char_th in range(1, len_seqs):
        for state in range(num_states):
            probability = omega[char_th - 1] + np.log(T[:, state]) + np.log(E[state, V[char_th]])
            probable_states[char_th - 1, state] = np.argmax(probability)
            omega[char_th, state] = np.max(probability)

    # Most observable sequences
    S = np.zeros(len_seqs)

    # Find the most probable last hidden state
    last_state = int(np.argmax(omega[-1, :]))
    S[-1] = last_state

    # Back tracking
    for i in range(2, len_seqs + 1):
        S[-i] = probable_states[len_seqs - 1 - i, last_state]
        last_state = int(probable_states[len_seqs - 1 - i, last_state])

    return S


'''
Build the occurrence matrix
'''


def build_occurrence_matrix(firstnames, surnames, seq, state_th_2_firstname_th, state_th_2_surname_th):

    M = np.zeros((len(firstnames), len(surnames)), int)
    probable_firstnames = []
    probable_surnames = []

    for s in seq:
        if s in state_th_2_firstname_th:
            probable_firstnames.append(state_th_2_firstname_th[s])
        elif s in state_th_2_surname_th:
            probable_surnames.append(state_th_2_surname_th[s])

    for i, j in zip(probable_firstnames, probable_surnames):
        M[i, j] += 1

    return M


if __name__ == '__main__':
    main()
