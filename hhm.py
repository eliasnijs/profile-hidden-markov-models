"""
Constructs a profile HMM given an input multiple sequence alignment, an alphabet, a pseudocount value and treshold. The function returns both the transition probabilities and emission probabilities.

This function uses a lot of branches and local variables. This is the case because I don't agree with the philosophy of splitting everything up in a lot of functions. Instead I try to write one, more coherent function. This, in my opinion, improves legibility as the user can read the function in a more linear way rather than jumping between a lot of smaller functions.

Parameters:
- treshold (float): the threshold used to determine which indices in the alignment correspond to insertions.
- pseudocount (float): the value added to each count before normalization.
- alphabet (str): a string representing the alphabet used in the sequences.
- alignment (list[str]): a list of strings representing the multiple sequence alignment.

Returns:
- A tuple of two dictionaries:
  - The first dictionary represents the transition probabilities of the profile HMM.
  - The second dictionary represents the emission probabilities of the profile HMM.

Example usage:
transition_probs, emission_probs = profile_HMM_pseudocounts(0.25, 0.04, 'ABCDE', ['ABD-DDB', 'CAACCEC'])
"""
def profile_HMM_pseudocounts(treshold: float, pseudocount: float, alphabet: str, alignment: list[str]):

    # Check which indices are insertions and determine the seed length
    insertion_indices = [j for j in range(len(alignment[0])) if sum(1 for i in range(len(alignment)) if alignment[i][j] == '-') / len(alignment) > treshold]
    seed_length = len(alignment[0]) - len(insertion_indices)

    # Construct profile counts
    emission_dict = {
        **{f"M{j}": {l: 0.0 for l in alphabet} for j in range(1, seed_length + 1)},
        **{f"I{j}": {l: 0.0 for l in alphabet} for j in range(0, seed_length + 1)},
    }

    transition_dict = {
        "S": {"I0": 0.0, "M1": 0.0, "D1": 0.0},
        "I0": {"I0": 0.0, "M1": 0.0, "D1": 0.0},
        **{f"M{j}": {f"I{j}": 0.0, f"M{j + 1}": 0.0, f"D{j + 1}": 0.0} for j in range(1, seed_length)},
        **{f"D{j}": {f"I{j}": 0.0, f"M{j + 1}": 0.0, f"D{j + 1}": 0.0} for j in range(1, seed_length)},
        **{f"I{j}": {f"I{j}": 0.0, f"M{j + 1}": 0.0, f"D{j + 1}": 0.0} for j in range(1, seed_length)},
        f"M{seed_length}": {f"I{seed_length}": 0.0, "E": 0.0},
        f"D{seed_length}": {f"I{seed_length}": 0.0, "E": 0.0},
        f"I{seed_length}": {f"I{seed_length}": 0.0, "E": 0.0},
    }

    for sequence in alignment:
        current_position = 'S'
        j_seed = 0
        for j, letter in enumerate(sequence):
            if j in insertion_indices:
                next_position = f"I{j_seed}"
                if letter != '-':
                    emission_dict[next_position][letter] += 1
                    transition_dict[current_position][next_position] += 1
                    current_position = next_position
            elif letter == '-':
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
    for emission_probs in emission_dict.values():
        total = sum(emission_probs.values())
        for letter, v in emission_probs.items():
            if total > 0:
                emission_probs[letter] = (v / total + pseudocount) / (1 + pseudocount * len(alphabet))
            else:
                emission_probs[letter] = pseudocount / (pseudocount * len(alphabet))

    for transition_probs in transition_dict.values():
        total = sum(transition_probs.values())
        for n2, v in transition_probs.items():
            if total > 0:
                transition_probs[n2] = (v / total + pseudocount) / (1 + pseudocount * len(transition_probs))
            else:
                transition_probs[n2] = pseudocount / (pseudocount * len(transition_probs))

    # Return the dictionaries
    return transition_dict, emission_dict
