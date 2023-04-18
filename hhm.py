import numpy as np
np.set_printoptions(linewidth=np.inf)


DELETION_CHARACHTER = '-'

def transition_dictionary_increment(dictionary, k1, k2):
        if k1 not in dictionary:
                dictionary[k1] = {}
        if k2 not in dictionary[k1]:
                dictionary[k1][k2] = 0
        dictionary[k1][k2] += 1

def emission_dictionary_increment(dictionary, k1, k2, alpha):
        if k1 not in dictionary:
                dictionary[k1] = {l: 0 for l in alpha}
        dictionary[k1][k2] += 1

def profile_HMM(treshold: float, pseudocount: float, alphabet: str, alignment: list[str]):
        # Construct dictionary mapping a letter onto an index for later use in our matrix
        alpha_to_index = {l:i for i,l in enumerate(alphabet)}

        # Construct the seed alignment and save the insertion indices
        insertion_indices = []
        for j in range(len(alignment[0])):
                spaces = 0
                for i in range(len(alignment)):
                        if alignment[i][j] == DELETION_CHARACHTER:
                                spaces += 1
                spaces_normalized = spaces / len(alignment)
                if (spaces_normalized > treshold):
                        insertion_indices.append(j)

        # Construct profile counts
        emission_dict = {}
        transition_dict = {}

        for i in range(len(alignment)):
                current_position = 'S'
                j_seed = 0
                for j, l in enumerate(alignment[i]):
                        next_position = ''
                        if l == DELETION_CHARACHTER:
                                next_position = f"D{j_seed+1}"
                                j_seed += 1
                        elif j in insertion_indices:
                                next_position = f"I{j_seed}"
                                emission_dictionary_increment(emission_dict, next_position, l, alphabet)
                        else:
                                next_position = f"M{j_seed+1}"
                                emission_dictionary_increment(emission_dict, next_position, l, alphabet)
                                j_seed += 1
                        transition_dictionary_increment(transition_dict, current_position, next_position)
                        current_position = next_position

        # Normalize profile counts
        for i, node in enumerate(emission_dict):
                totalcount = 0
                for letter in emission_dict[node]:
                        totalcount += emission_dict[node][letter]
                for letter in emission_dict[node]:
                        emission_dict[node][letter] /= totalcount

        for i, node1 in enumerate(transition_dict):
                totalcount = 0
                for node2 in transition_dict[node1]:
                        totalcount += transition_dict[node1][node2]
                for node2 in transition_dict[node1]:
                        transition_dict[node1][node2] /= totalcount

        print(emission_dict)
        print(transition_dict)


def profile_HMM_2(treshold: float, pseudocount: float, alphabet: str, alignment: list[str]):
        # Construct dictionary mapping a letter onto an index for later use in our matrix
        alpha_to_index = {l:i for i,l in enumerate(alphabet)}

        # Construct the seed alignment and save the insertion indices
        insertion_indices = []
        for j in range(len(alignment[0])):
                spaces = 0
                for i in range(len(alignment)):
                        if alignment[i][j] == DELETION_CHARACHTER:
                                spaces += 1
                spaces_normalized = spaces / len(alignment)
                if (spaces_normalized > treshold):
                        insertion_indices.append(j)
        seed_length = len(alignment[0]) - len(insertion_indices)

        # Construct profile counts
        transition_matrix = np.zeros((3, 3 * seed_length + 2))
        emission_matrix = np.zeros((len(alphabet), 2 * seed_length + 1))
        string = alignment[0]

        for i in range(len(alignment)):
                print(alignment[i])
                row = 1
                j_seed = 0
                for j, l in enumerate(alignment[i]):
                        new_row = 0
                        column = j_seed * 3
                        if j in insertion_indices:
                                if l != DELETION_CHARACHTER:
                                        emission_matrix[alpha_to_index[l]][j_seed + seed_length] += 1
                                new_row = 2
                        elif l == DELETION_CHARACHTER:
                                new_row = 1
                                column += 2
                                j_seed += 1
                        else:
                                emission_matrix[alpha_to_index[l]][j_seed] += 1
                                new_row = 0
                                column += 1
                                j_seed += 1
                        transition_matrix[row, column] += 1
                        row = new_row
                j = len(alignment[i])
                while j in insertion_indices:
                        transition_matrix[row, j_seed * 3] += 1
                        row = 2
                transition_matrix[row, j_seed * 3 + 1] += 1

        transition_matrix += pseudocount
        emission_matrix += pseudocount

        print(emission_matrix)

        # Normalize profile counts
        for j in range(0, transition_matrix.shape[1], 3):
                for i in range(3):
                        transition_matrix[i][j:j+3] /= np.sum(transition_matrix[i][j:j+3])
        for j in range(0, emission_matrix.shape[1]):
                emission_matrix[:,j] /= np.sum(emission_matrix[:,j])

        print(transition_matrix)
        print(emission_matrix)
        emission_dict = 0
        return transition_matrix, emission_dict


# theta = 0.35
# pseudocount = 0.04
# alphabet = 'ACDEF'
# alignment = [ 'ACDEFACADF', 'AFDA---CCF', 'A--EFD-FDC', 'ACAEF--A-C', 'ADDEFAAADF']
theta = 0.25
pseudocount = 0.04
alphabet = 'ABCDE'
alignment = ['ABD-DDB', 'CAACCEC']

profile_HMM_2(theta, pseudocount, alphabet, alignment)
