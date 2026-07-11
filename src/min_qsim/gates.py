import numpy as np

"""Single-Qubit Gates"""
#Constant Gates
#-------------------------------------------------
I = np.array([[1, 0],
              [0, 1]])

X = np.array([[0, 1], 
              [1, 0]])

Y = np.array([[0, -1j], 
              [1j, 0]])

Z = np.array([[1, 0],
              [0, -1]])

H = np.array([[1, 1],
              [1, -1]]) / np.sqrt(2)

S = np.array([[1, 0],
              [0, 1j]])

T = np.array([[1, 0],
              [0, np.exp(1j * np.pi / 4)]])

CNOT = np.array([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1],
                 [0, 0, 1, 0]])

CZ = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, -1]])

SWAP = np.array([[1, 0, 0, 0],
                 [0, 0, 1, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1]])

CCX = np.array([[1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 0, 1, 0]])
#-------------------------------------------------

#Parameterised Gates
#-------------------------------------------------
def rx(theta):
    """
    Intializes a rotationary matrix about the x-axis with user given angle

    Args:
        theta: angle to turn

    Return:
        Rotionary matrix about the x-axis with angle theta
    """
    return np.array([[np.cos(theta/2), -1j * np.sin(theta/2)],
                     [-1j * np.sin(theta/2), np.cos(theta/2)]])

def ry(theta):
    """
    Intializes a rotationary matrix about the y-axis with user given angle

    Args:
        theta: angle to turn

    Return: 
        Rotationary matrix about the y-axis with angle theta
    """
    return np.array([[np.cos(theta/2), -np.sin(theta/2)],
                     [np.sin(theta/2), np.cos(theta/2)]])

def rz(theta):
    """
    Intializes a rotationary matrix about the z-axis with user given angle

    Args:
        theta: angle to turn

    Return:
        Rotationary matrix about the z-axis with angle theta
    """
    return np.array([[np.exp(-1j * theta/2), 0],
                     [0, np.exp(1j * theta/2)]])

def apply_gate(state, gate, qubit_index: int, n: int):
    """
    Applies a single-qubit gate to the qubit at qubit_index in an n-qubit system.

    Args:
        state: numpy array, length 2**n, the current statevector
        gate: 2x2 numpy array, the single-qubit gate to apply
        qubit_index: which qubit (0-indexed) the gate acts on
        n: total number of qubits in the system

    Returns:
        new_state: numpy array, length 2**n, state after applying the gate

    Raises:
        ValueError: if gate is not 2x2
        ValueError: if qubit_index >= n or qubit_index < 0
        ValueError: if gate is not unitary
    """
    if gate.shape != (2, 2): raise ValueError("Gate is not 2x2")
    if qubit_index >= n or qubit_index < 0: raise ValueError("Invalid qubit_index")
    if not np.allclose(gate.conj().T @ gate, np.eye(2)): raise ValueError("Gate is not unitary")

    #Generating gate
    full_matrix = np.array([[1]])
    for i in range(n):
        if i != qubit_index:
            full_matrix = np.kron(full_matrix, I)
        else:
            full_matrix = np.kron(full_matrix, gate)
    
    new_state = full_matrix @ state
    return new_state
#-------------------------------------------------

"""Multi-Qubit Gates"""

def apply_cnot(state, control: int, target: int, n: int):
    """
    Applies a CNOT gate to a qubit

    Args:
        state: n-qubit system
        control: the contrl qubit 
        target: the target qubit
        n: total number of qubits in the system

    Returns:
        new_state: numpy array, length 2**n, state after applying CNOT to target qubit
    
    Raises:
        ValueError: if control >= n or control < 0
        ValueError: if state is not the right length
        ValueError: if target == control 
        ValueError: if target >= n or target < 0
    """

    if control >= n or control < 0: raise ValueError("Control out of range")
    if target >= n or target < 0: raise ValueError("Target out of range")
    if target == control: raise ValueError("Target and control cannot be the same")
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")

    control_offset = n - control - 1
    target_offset = n - target - 1

    new_state = np.zeros(2 ** n, dtype=complex)
    for i in range(2 ** n):
        control_bit = (i >> control_offset) & 1
        if control_bit == 0:
            new_state[i] = state[i]
        else:
            i_flipped = i ^ (1 << target_offset)
            new_state[i_flipped] = state[i]

    return new_state

def apply_cz(state, qubit_a: int, qubit_b: int, n: int):
    """
    Applies a CZ gate to a qubit

    Args:
        state: n-qubit system
        qubit_a: the qubit_a qubit (Qubit A and B are symmetric, swapping them yields the same results)
        qubit_b: the qubit_b qubit
        n: total number of qubits in the system
    
    Returns:
        new_state: numpy array, length 2**n, state after applying CZ to target qubit
    
    Raises:
        ValueError: if qubit_a >= n or qubit_a < 0
        ValueError: if qubit_b >= n  or qubit_b < 0
        ValueError: if qubit_b == qubit_a
        ValueError: if state is not the right length
    """

    if qubit_a >= n or qubit_a < 0: raise ValueError("Qubit A out of range")
    if qubit_b >= n or qubit_b < 0: raise ValueError("Qubit B out of range")
    if qubit_b == qubit_a: raise ValueError("Qubit A and B cannot be the same")
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")

    a_offset = n - qubit_a - 1
    b_offset = n - qubit_b - 1

    new_state = np.zeros(2 ** n, dtype=complex)
    for i in range(2 ** n):
        a_bit = (i >> a_offset) & 1
        b_bit = (i >> b_offset) & 1
        if a_bit == b_bit == 1:
            new_state[i] = -state[i]
        else:
            new_state[i] = state[i]
    
    return new_state

def apply_swap(state, qubit_a: int, qubit_b: int, n: int):
    """
    Applies a SWAP gate to a qubit

    Args:
        state: n-qubit system
        qubit_a: swaps with qubit_b (qubit_a and qubit_b are the symmetric, changing them yields the same result)
        qubit_b: swaps with qubit_a
        n: total number of qubits in the system

    Returns:
        new_state: numpy array, length 2**n, state after applying SWAP

    Raises:
        ValueError: if qubit_a >= n or qubit_a < 0
        ValueError: if qubit_b >= n  or qubit_b < 0
        ValueError: if qubit_b == qubit_a
        ValueError: if state is not the right length
    """
    if qubit_a >= n or qubit_a < 0: raise ValueError("Qubit A out of range")
    if qubit_b >= n or qubit_b < 0: raise ValueError("Qubit B out of range")
    if qubit_a == qubit_b: raise ValueError("Qubit A and Qubit B cannot be the same")
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")

    a_offset = n - qubit_a - 1
    b_offset = n - qubit_b - 1

    new_state = np.zeros(2 ** n, dtype=complex)
    for i in range(2 ** n):
        a_bit = (i >> a_offset) & 1
        b_bit = (i >> b_offset) & 1
        if a_bit == b_bit:
            new_state[i] = state[i]
        else:
            combined_mask = (1 << a_offset) | (1 << b_offset)
            ab_flipped = i ^ combined_mask
            new_state[ab_flipped] = state[i]

    return new_state

def apply_ccx(state, qubit_a: int, qubit_b: int, target: int, n: int):
    """
    Applies a CCX gate to a qubit at the given control location

    Args:
        state: n-qubit system
        qubit_a: input into CCX gate (qubit_a and qubit_b are symmetric, swapping them yield the same result)
        qubit_b: input into CCX gate
        target: flips if qubit_a == qubit_b == 1
        n: total number of qubits in the system

    Returns:
        new_state: numpy array, length 2**n, state after applying CCX

    Raises:
        ValueError: if qubit_a > n or qubit_a < 0
        ValueError: if qubit_b > n or qubit_b < 0
        ValueError: if target > n or target < 0
        ValueError: if qubit_a == qubit_b or qubit_a == target or qubit_b == target
        ValueError: if state is not the right length
        ValueError: if state < 2**3
        ValueError: if n < 3
    """
    if qubit_a >= n or qubit_a < 0: raise ValueError("Qubit A out of range")
    if qubit_b >= n or qubit_b < 0: raise ValueError("Qubit B out of range")
    if target >= n or target < 0: raise ValueError("Target out of range")
    if qubit_a == target or qubit_a == qubit_b or qubit_b == target: raise ValueError("Qubit A, B and target must be unique")
    if n < 3 or len(state) < 2**3: raise ValueError("State vector too small")
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")

    a_offset = n - qubit_a - 1
    b_offset = n - qubit_b - 1
    target_offset = n - target - 1

    new_state = np.zeros(2 ** n, dtype=complex)
    for i in range(2 ** n):
        a_bit = (i >> a_offset) & 1
        b_bit = (i >> b_offset) & 1
        if a_bit == b_bit == 1:
            target_flipped = i ^ (1 << target_offset)
            new_state[target_flipped] = state[i]
        else:
            new_state[i] = state[i]
    
    return new_state