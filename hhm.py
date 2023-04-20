import numpy as np

np.set_printoptions(linewidth=np.inf)

DELETION_CHARACHTER = '-'


def profile_HMM_pseudocounts(treshold: float, pseudocount: float, alphabet: str, alignment: list[str]):
    # Construct the seed alignment and save the insertion indices
    insertion_indices = [j for j in range(len(alignment[0])) if sum(1 for i in range(len(alignment)) if alignment[i][j] == DELETION_CHARACHTER) / len(alignment) > treshold]
    seed_length = len(alignment[0]) - len(insertion_indices)

    # Construct profile counts
    emission_dict = {"I0": {l: 0.0 for l in alphabet}}
    for j in range(1, seed_length + 1):
        emission_dict[f"M{j}"] = {l: 0.0 for l in alphabet}
        emission_dict[f"I{j}"] = {l: 0.0 for l in alphabet}

    transition_dict = {
        "S": {"I0": 0.0, "M1": 0.0, "D1": 0.0},
        "I0": {"I0": 0.0, "M1": 0.0, "D1": 0.0},
        f"M{seed_length}": {f"I{seed_length}": 0.0, "E": 0.0},
        f"D{seed_length}": {f"I{seed_length}": 0.0, "E": 0.0},
        f"I{seed_length}": {f"I{seed_length}": 0.0, "E": 0.0}
    }
    for j in range(1, seed_length):
        transition_dict[f"M{j}"] = {f"I{j}": 0.0, f"M{j + 1}": 0.0, f"D{j + 1}": 0.0}
        transition_dict[f"D{j}"] = {f"I{j}": 0.0, f"M{j + 1}": 0.0, f"D{j + 1}": 0.0}
        transition_dict[f"I{j}"] = {f"I{j}": 0.0, f"M{j + 1}": 0.0, f"D{j + 1}": 0.0}

    for sequence in alignment:
        current_position = 'S'
        j_seed = 0
        for j, letter in enumerate(sequence):
            if j in insertion_indices:
                next_position = f"I{j_seed}"
                if letter != DELETION_CHARACHTER:
                    emission_dict[next_position][letter] += 1
                    transition_dict[current_position][next_position] += 1
                    current_position = next_position
            elif letter == DELETION_CHARACHTER:
                next_position = f"D{j_seed + 1}"
                transition_dict[current_position][next_position] += 1
                current_position = next_position
                j_seed += 1
            else:
                next_position = f"M{j_seed + 1}"
                emission_dict[next_position][letter] += 1
                transition_dict[current_position][next_position] += 1
                current_position = next_position
                j_seed += 1
        j = len(sequence)
        while j in insertion_indices:
            transition_dict[current_position][f"I{j_seed}"] += 1
            current_position = f"I{j_seed}"
        transition_dict[current_position]["E"] += 1

    # Normalize profile counts and factor in pseudocounts
    for i, n in enumerate(emission_dict):
        total = sum(emission_dict[n].values())
        for letter, v in emission_dict[n].items():
            if total > 0:
                emission_dict[n][letter] = (v / total + pseudocount) / (1 + pseudocount * len(alphabet))
            else:
                emission_dict[n][letter] = pseudocount / (pseudocount * len(alphabet))

    for n1 in transition_dict:
        total = sum(transition_dict[n1].values())
        for n2, v in transition_dict[n1].items():
            if total > 0:
                transition_dict[n1][n2] = (v / total + pseudocount) / (1 + pseudocount * len(transition_dict[n1]))
            else:
                transition_dict[n1][n2] = pseudocount / (pseudocount * len(transition_dict[n1]))

    # Return the dictionaries
    return transition_dict, emission_dict


# theta = 0.35
# pseudocount = 0.04
# alphabet = 'ACDEF'
# alignment = [ 'ACDEFACADF', 'AFDA---CCF', 'A--EFD-FDC', 'ACAEF--A-C', 'ADDEFAAADF']
theta = 0.25
pseudocount = 0.04
alphabet = 'ABCDE'
alignment = ['ABD-DDB', 'CAACCEC']

profile_HMM_pseudocounts(theta, pseudocount, alphabet, alignment)
