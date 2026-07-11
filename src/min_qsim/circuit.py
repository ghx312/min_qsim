from min_qsim import state, gates, measurement, debug
import numpy as np

class Circuit:
    def __init__(self, n: int, tol: float=1e-6):
        """
        Initializes an n-qubit circuit in the |00...0> state.

        Args:
            n: number of qubits
            state: statevector
            tol: tolerance

        """
        self.n = n
        self.state = state.init_state(n)
        self.tol = tol

    def set_state(self, vector):
        """
        Replaces the circuit's current statevector with a user-provided one.

        Args:
            vector: array-like, length 2**n, the new statevector to load

        Returns:
            self, to allow method chaining
        """
        self.state = state.custom_state(vector, self.n)
        return self

    def x(self, qubit_index: int):
        """
        Applies a X gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply X to

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.X, qubit_index, self.n)
        return self
    
    def y(self, qubit_index: int):
        """
        Applies a Y gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply Y to

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.Y, qubit_index, self.n)
        return self
    
    def z(self, qubit_index: int):
        """
        Applies a Z gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply Z to

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.Z, qubit_index, self.n)
        return self
    
    def s(self, qubit_index: int):
        """
        Applies a S gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply S to

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.S, qubit_index, self.n)
        return self

    def t(self, qubit_index: int):
        """
        Applies a T gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply T to

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.T, qubit_index, self.n)
        return self

    def h(self, qubit_index: int):
        """
        Applies a Hadamard gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply H to

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.H, qubit_index, self.n)
        return self
    
    def rx(self, qubit_index: int, theta: float):
        """
        Applies a RX gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply RX to
            theta: the angle to turn about the x-axis

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.rx(theta), qubit_index, self.n)
        return self
    
    def ry(self, qubit_index: int, theta: float):
        """
        Applies a RY gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply RY to
            theta: the angle to turn about the y-axis

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.ry(theta), qubit_index, self.n)
        return self
    
    def rz(self, qubit_index: int, theta: float):
        """
        Applies a RZ gate to the specified qubit.

        Args:
            qubit_index: the qubit to apply RZ to
            theta: the angle to turn about the z-axis

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_gate(self.state, gates.rz(theta), qubit_index, self.n)
        return self

    def cnot(self, control, target):
        """
        Applies a CNOT gate according to control to target

        Args:
            control: control qubit index
            target: target qubit index

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_cnot(self.state, control, target, self.n)
        return self
    
    def cz(self, qubit_a: int, qubit_b: int):
        """
        Applies a CZ gate according to qubit_a and qubit_b

        Args:
            qubit_a: the qubit_a qubit (Qubit A and B are symmetric, swapping them yields the same results)
            qubit_b: the qubit_b qubit

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_cz(self.state, qubit_a, qubit_b, self.n)
        return self
    
    def swap(self, qubit_a: int, qubit_b: int):
        """
        Applies a SWAP gate according to qubit_a and qubit_b

        Args:
            qubit_a: the qubit_a qubit (Qubit A and B are symmetric, swapping them yields the same results)
            qubit_b: the qubit_b qubit

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_swap(self.state, qubit_a, qubit_b, self.n)
        return self
    
    def ccx(self, qubit_a: int, qubit_b: int, target: int):
        """
        Applies a CCX gate according to qubit_a and qubit_b to target

        Args:
            qubit_a: the qubit_a qubit (Qubit A and B are symmetric, swapping them yields the same results)
            qubit_b: the qubit_b qubit
            target: target qubit index

        Returns:
            self, to allow method chaining
        """
        self.state = gates.apply_ccx(self.state, qubit_a, qubit_b, target, self.n)
        return self
    
    def fmeasure(self):
        """
        Performs a full measurement, collapsing the entire state.

        Returns:
            observed: the basis string outcome
        """
        self.state, new_basis = measurement.full_measurement(self.state, self.n)
        return new_basis

    def pmeasure(self, qubits: list):
        """
        Performs a partial measurement on the specified qubits.

        Args:
            qubits: list of qubit indices to measure

        Returns:
            observed: the combo string outcome for the measured qubits
        """
        self.state, measured_basis = measurement.partial_measurement(self.state, self.n, qubits)
        return measured_basis
    
    def reset(self, qubit_index: int):
        """
        Resets the specified qubit to |0>, collapsing any entanglement.

        Args:
            qubit_index: resetted qubits

        Returns:
            self, to allow method chaining
        """
        self.state = measurement.reset_qubit(self.state, self.n, qubit_index)
        return self

    def probabilities(self, basis_state: str = None):
        """
        Args:
            basis_state (Optional): only possibility of basis state
        Returns the probabilities of basis states without collapsing self.state.
        """
        prob = measurement.get_probabilities(self.state, self.n, basis_state)
        return prob
    
    def compose(self, other):
        """
        Combines this circuit with another independent circuit via tensor
        product, producing a combined state over both circuits' qubits.

        Args:
            other: another Circuit instance

        Returns:
            self, with self.n and self.state updated to reflect the combined system
            (self.n becomes self.n + other.n; self.state becomes the kron product)
        """
        self.state = np.kron(self.state, other.state)
        self.n = self.n + other.n
        return self

    def __repr__(self):
        """
        Args:
            tol: tolerance set to 1e-6 default
        
        Returns a readable string representation of the circuit's current state.
        """
        state_string = debug.state_to_string(self.state, self.n, self.tol)
        return state_string