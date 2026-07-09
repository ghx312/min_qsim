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

def apply_gate(state, gate, qubit_index, n):
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

def apply_cnot(state, control, n):
    """
    Applies a CNOT gate to a qubit at the given control location

    Args:
        state: n-qubit system
        control: the contrl qubit indicates if the next qubit will flip
        n: total number of qubits in the system

    Returns:
        new_state: numpy array, length 2**n, state after applying CNOT at control
    
    Raises:
        ValueError: if control + 1 >= n or control < 0
        ValueError: if state is not the right length
    """

    if control + 1 >= n or control < 0: raise ValueError("Control out of range")
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")

    full_matrix = np.array([[1]])
    for i in range(n):
        if i == control: full_matrix = np.kron(full_matrix, CNOT)
        elif i == control + 1: continue
        else: full_matrix = np.kron(full_matrix, I)
    
    new_state = full_matrix @ state
    return new_state

def apply_cz(state, control, n):
    """
    Applies a CZ gate to a qubit at the given control location

    Args:
        state: n-qubit system
        control: the control qubit and the control + 1 qubit determines the changes made by CZ gate
        n: total number of qubits in the system
    
    Returns:
        new_state: numpy array, length 2**n, state after applying CZ at control
    
    Raises:
        ValueError: if control + 1 >= n or control < 0
        ValueError: if state is not the right length
    """

    if control + 1 >= n or control < 0: raise ValueError("Control out of range")
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")

    full_matrix = np.array([[1]])
    for i in range(n):
        if i == control: full_matrix = np.kron(full_matrix, CZ)
        elif i == control + 1: continue
        else: full_matrix = np.kron(full_matrix, I)
    
    new_state = full_matrix @ state
    return new_state

def apply_swap(state, control, n):
    """
    Applies a SWAP gate to a qubit at the given control location

    Args:
        state: n-qubit system
        control: swaps the control and control + 1 qubit
        n: total number of qubits in the system

    Returns:
        new_state: numpy array, length 2**n, state after applying SWAP at control

    Raises:
        ValueError: if control + 1 >= n or control < 0
        ValueError: if state is not the right length
    """
    if control + 1 >= n or control < 0: raise ValueError("Control out of range")
    if len(state) != 2 ** n: raise ValueError("Incorrect values provided")

    full_matrix = np.array([[1]])
    for i in range(n):
        if i == control: full_matrix = np.kron(full_matrix, SWAP)
        elif i == control + 1: continue
        else: full_matrix = np.kron(full_matrix, I)
    
    new_state = full_matrix @ state
    return new_state

def apply_toffoli(state, control, n):
    """
    Applies a CCX/toffoli gate to a qubit at the given control location

    Args:
        state: n-qubit system
        control: swaps the control and control + 1 qubit
        n: total number of qubits in the system

    Returns:
        new_state: numpy array, length 2**n, state after applying SWAP at control

    Raises:
        ValueError: if control + 1 >= n or control < 0
        ValueError: if state is not the right length
    """
    pass
