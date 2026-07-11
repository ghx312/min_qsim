import numpy as np

def init_state(n: int): #n is the number of qubits in the system
    """
    Initialize an n-qubit system to the |00...0> state.
    
    Args: 
        n: number of qubits

    Returns:
        A normalized numpy array of length 2**n, with 1 at index 0 and 0 elsewhere.
    """
    if n <= 0: raise ValueError("n must be at least 1")
    
    vector = np.zeros(2 ** n, dtype=complex)
    vector[0] = 1
    return vector

def custom_state(vector, n: int):
    """
    Accept a user-provided vector as a custum quantum state

    Args:
        vector: array-like of complex/float amplitudes [List]
        n: total number of qubits in the system

    Returns:
        A numpy array representing a validated state.

    Raises:
        ValueError: if len(vector) != 2 ** n
        ValueError: if the vector is not normalized
    """
    vector = np.array(vector)

    vector_length = len(vector)
    if vector_length == 0: raise ValueError("Invalid vector length")
    if vector_length != 2 ** n: raise ValueError("Incorrect values provided")
    
    squared_magnitude = np.abs(vector) ** 2
    total_probability = np.sum(squared_magnitude)
    is_normalized = np.isclose(total_probability, 1.0)
    if not is_normalized: raise ValueError("Statevector not normalized")

    return vector