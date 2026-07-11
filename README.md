# min_qsim

A minimalistic n-qubit statevector quantum circuit simulator, built from scratch in NumPy for educational purposes. No external quantum computing libraries — every gate, measurement, and state operation is implemented directly on top of the raw statevector using explicit bit-manipulation, not black-box matrix libraries.

This project exists to make the bridge between pen-and-paper quantum mechanics (bra-ket notation, unitary matrices, tensor products, Born's rule) and working code as transparent as possible.

## Installation

```bash
pip install min-qsim
```

Requires Python >= 3.9 and NumPy >= 1.24.

## Quick example

```python
import min_qsim as q

# Bell state: |00> -> H(q0) -> CNOT(q0, q1)
qc = q.Circuit(2)
qc.h(0).cnot(0, 1)
print(qc)
# 0.707|00> + 0.707|11>

outcome = qc.fmeasure()
print(outcome)
# '00' or '11' -- collapses the entangled pair together
```

## The `Circuit` class

`Circuit` wraps every operation in the library so you don't have to manually thread a statevector and qubit count through function calls. Every gate method mutates the circuit's internal state and returns `self`, so calls can be chained.

```python
qc = q.Circuit(n, tol=1e-6)  # n qubits, initialized to |00...0>
```

- `tol` controls how small an amplitude has to be before it's hidden from the printed representation (`print(qc)`), not the internal precision of any calculation.

### Single-qubit gates
```python
qc.x(qubit)          # Pauli-X
qc.y(qubit)          # Pauli-Y
qc.z(qubit)          # Pauli-Z
qc.h(qubit)          # Hadamard
qc.s(qubit)          # S (phase) gate
qc.t(qubit)          # T gate
qc.rx(qubit, theta)  # rotation about x-axis
qc.ry(qubit, theta)  # rotation about y-axis
qc.rz(qubit, theta)  # rotation about z-axis
```

### Multi-qubit gates
```python
qc.cnot(control, target)
qc.cz(qubit_a, qubit_b)          # symmetric: order doesn't matter
qc.swap(qubit_a, qubit_b)        # symmetric: order doesn't matter
qc.ccx(qubit_a, qubit_b, target) # Toffoli/CCX; qubit_a/qubit_b symmetric
```

All multi-qubit gates work on **arbitrary, non-adjacent qubits** — there's no requirement that control and target be next to each other.

### Measurement
```python
qc.fmeasure()          # full measurement: collapses the entire state, returns a basis string like '011'
qc.pmeasure([0, 2])    # partial measurement on specific qubits, returns a combo string
qc.probabilities()     # dict of {basis_state: probability} without collapsing anything
qc.probabilities('01') # probability of one specific basis state
```

### Other operations
```python
qc.reset(qubit)              # forces a qubit back to |0>, collapsing any entanglement first
qc.set_state(vector)          # load a custom statevector (must match 2**n in length and be normalized)
qc.compose(other_circuit)     # tensor-product this circuit with another; self's qubits come first
```

### Inspecting a circuit
```python
print(qc)   # human-readable ket notation, e.g. "0.707|00> + 0.707|11>"
```

## Function-level API

Every `Circuit` method is a thin wrapper around standalone functions, all of which are also importable directly if you'd rather work with raw statevectors:

```python
from min_qsim import (
    init_state, custom_state,
    X, Y, Z, H, S, T, CNOT, CZ, SWAP, CCX,   # gate constants
    rx, ry, rz,                               # parameterized gate constructors
    apply_gate, apply_cnot, apply_cz, apply_swap, apply_ccx,
    get_probabilities, sampling, full_measurement, partial_measurement,
    reset_qubit, is_normalized,
    state_to_string, format_amplitude,
)
```

- `init_state(n)` / `custom_state(vector, n)` — build a statevector from scratch or from a user-provided array.
- `apply_gate(state, gate, qubit_index, n)` — apply any single-qubit gate to any qubit in an n-qubit system.
- `apply_cnot` / `apply_cz` / `apply_swap` / `apply_ccx` — multi-qubit gates, all generalized to arbitrary qubit positions via bit-manipulation on basis-state indices.
- `get_probabilities(state, n, basis_state=None)` — raw probabilities without collapsing the state.
- `sampling(state, n, shots)` — repeated measurement, returns frequency counts.
- `full_measurement` / `partial_measurement` — collapse and observe (single draw, not multi-shot).
- `reset_qubit(state, n, qubit_index)` — measure-and-correct a single qubit back to |0>.
- `is_normalized(state, n)` — sanity-check a statevector sums to probability 1.
- `state_to_string(state, n, tol=1e-6)` — human-readable ket notation for any statevector.

## Design notes

- **Leftmost-first convention**: qubit 0 is the most-significant bit of every basis-state index.
- Multi-qubit gates use direct bit-manipulation (`(i >> (n-1-q)) & 1` to extract a bit, XOR/OR of masks to flip) rather than constructing full kron-chain matrices — this keeps memory and computation manageable and mirrors how you'd reason about the gate by hand.
- `partial_measurement`/`reset_qubit` preserve relative phase and magnitude between surviving amplitudes; only a fully-collapsed single-basis-state result can have its amplitude hardcoded to 1.0 (global phase is unobservable).

## License

MIT
