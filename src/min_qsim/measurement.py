import numpy as np

def get_probabilities(state, n: int, basis_state: str=None):
    """
    Get the probabilities of the different basis_state in your statevector, pass in basis_state
    if they want only 1 probability

    Args:
        state: n-qubit system
        n: total number of qubits in the system
        basis_state: string of basis_state if they want only this state's probability

    Returns:
        if basis_state == None: returns a dictionary containing the proability of every basis state
        if basis_state != None: returns the probability of that state as a float

    Raises:
        ValueError: state is not right length
        ValueError: statevector is not normalized
        ValueError: basis_state does not exist
    """
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")
    squared_magnitude = np.abs(state) ** 2
    total_probability = np.sum(squared_magnitude)
    is_normalized = np.isclose(total_probability, 1.0)
    if not is_normalized: raise ValueError("Statevector not normalized")

    prob_dict = {}
    for i in range(2 ** n):
        current_state = bin(i)[2:].zfill(n)
        current_prob = round(float(squared_magnitude[i]), 3)
        prob_dict.update({current_state : current_prob})
    
    if basis_state == None:
        return prob_dict    
    elif basis_state in prob_dict.keys():
        return prob_dict[basis_state]
    else:
        raise ValueError("Invalide Basis State")
    
def sampling(state, n: int, shots: int):
    """
    Measures and records the results of a qubit system shot number of times

    Args:
        state: n-qubit system
        n: total number of qubits in the system
        shots: Number of times the measurement is ran and the results is collected

    Returns:
        A dictionary containing frequency of basis vectors

    Raises:
        ValueError: state is not right length
        ValueError: statevector is not normalized
    """
    prob_dict = get_probabilities(state, n)
    basis_states = list(prob_dict.keys())
    probability = list(prob_dict.values())
    data_dict = dict.fromkeys(basis_states, 0)
    for _ in range(shots):
        idx = np.random.choice(len(basis_states), p=probability)
        data_dict[basis_states[idx]] += 1

    return data_dict
        
def full_measurement(state, n: int):
    """
    Fully measures a qubit and outputs its measured basis

    Args:
        state: n-qubit system
        n: total number of qubits in the system
    
    Returns:
        Final collapsed state + Measured basis
    
    Raises:
    ValueError: state is not right length
    ValueError: statevector is not normalized
    """
    data = sampling(state, n, 1)
    freq_idx = 0
    freq_data = list(data.values())
    for i in range(len(freq_data)):
        if freq_data[i] == 1:
            freq_idx = i
    basis_states = list(data.keys())

    observed_basis = basis_states[freq_idx]
    observed_idx = int(observed_basis, 2)

    new_state = np.zeros(2 ** n, dtype=complex)
    new_state[observed_idx] = 1.0

    return new_state, observed_basis

def partial_measurement(state, n: int, qubit_to_measure):
    """
    Takes partial measurements and returns a normalized statevector of possible states
    after partial measurmeents
    
    Args:
        state: n-qubit system
        n: total number of qubits in the system
        qubit_to_measure: A list consisting of the index of qubits the user wants to measure

    Returns:
        A normalized statevector dependant on the qubit_dict provided and state

    Raises:
        ValueError: Initial statevector not normalized
        ValueError: len(state) != 2 ** n
        ValueError: Index in qubit_to_measure out of range
    """ 
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")
    squared_magnitude = np.abs(state) ** 2
    total_probability = np.sum(squared_magnitude)
    is_normalized = np.isclose(total_probability, 1.0)
    if not is_normalized: raise ValueError("Statevector not normalized")
    for i in range(len(qubit_to_measure)):
        if qubit_to_measure[i] >= n or qubit_to_measure[i] < 0:
            raise ValueError("Qubit index out of range")

    #Constructing qubits_to_measure statevector
    m_qubit_probabilities = np.zeros(2 ** len(qubit_to_measure), dtype=float)
    for i in range(2 ** n):
        current_state = bin(i)[2:].zfill(n)
        current_state_int = 0
        for j in range(len(qubit_to_measure)):
            add = current_state[qubit_to_measure[j]]
            if add == "1": current_state_int += 2 ** (len(qubit_to_measure) - j - 1)
        m_qubit_probabilities[current_state_int] += np.abs(state[i]) ** 2
    
    #Normalizing m_qubit_probabilities
    m_qubit_probabilities /= np.sum(m_qubit_probabilities)

    #Sampling over selected qubit's probabilities
    m_basis_state = [bin(i)[2:].zfill(len(qubit_to_measure)) for i in range(2 ** len(qubit_to_measure))]
    basis_state_idx = np.random.choice(len(m_basis_state), p=m_qubit_probabilities)
    measured = m_basis_state[basis_state_idx]
    targets = [int(bit) for bit in measured]

    #Constructing new normalized statevector given target values
    new_state = np.zeros(2 ** n, dtype=complex)
    for i in range(2 ** n):
        current_state = bin(i)[2:].zfill(n)
        match = True
        for q, t in zip(qubit_to_measure, targets):
                offset = n - q - 1
                if ((i >> offset) & 1) != t:
                    match = False
                    break
        if match: new_state[i] = state[i]

    new_state /= np.sqrt(m_qubit_probabilities[basis_state_idx])
    return new_state