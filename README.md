# min_qsim

A minimalistic n-qubit statevector quantum circuit simulator, built from scratch in NumPy for educational purposes. No external quantum computing libraries — every gate, state, and measurement routine is implemented directly on top of the raw statevector.

## Installation

```bash
git clone https://github.com/yourhandle/min_qsim.git
cd min_qsim
pip install -e .
```

Requires Python >= 3.9 and NumPy >= 1.24.

## Quick example

```python
import min_qsim as q

# Bell state: |00> -> H(q0) -> CNOT(q0, q1)
state = q.init_state(2)
state = q.apply_gate(state, q.H, 0, 2)
state = q.apply_cnot(state, 0, 1, 2)

print(q.get_probabilities(state, 2))
# {'00': 0.5, '01': 0.0, '10': 0.0, '11': 0.5}

print(q.sampling(state, 2, shots=1000))
# {'00': 502, '01': 0, '10': 0, '11': 498}
```

## Status

### Completed

**State initialization** (`state.py`)
- `init_state(n)` — initializes an n-qubit register to |00...0>
- `custom_state(vector)` — validates and loads an arbitrary user-provided statevector (checks power-of-2 length and normalization)

**Gates** (`gates.py`)
- Constant single-qubit gates: `I`, `X`, `Y`, `Z`, `H`, `S`, `T`
- Constant multi-qubit gates (as raw matrices): `CNOT`, `CZ`, `SWAP`, `CCX`
- Parameterised rotation gates: `rx(theta)`, `ry(theta)`, `rz(theta)`
- `apply_gate(state, gate, qubit_index, n)` — applies any single-qubit gate to a target qubit via full tensor-product expansion, with unitarity/shape/index validation
- `apply_cnot(state, control, target, n)` — bitmask-based CNOT application (no full matrix construction)
- `apply_cz(state, qubit_a, qubit_b, n)` — bitmask-based CZ application
- `apply_swap(state, qubit_a, qubit_b, n)` — bitmask-based SWAP application
- `apply_toffoli(state, qubit_a, qubit_b, target, n)` — bitmask-based CCX application

**Measurement** (`measurement.py`)
- `get_probabilities(state, n, basis_state=None)` — returns the probability of every basis state, or a single basis state's probability
- `sampling(state, n, shots)` — repeatedly samples the distribution to build up measurement statistics without collapsing state between shots
- `full_measurement(state, n)` — performs a single full projective measurement, returning the collapsed state and the observed basis string
- `partial_measurement(state, n, qubit_to_measure)` — measures a subset of qubits, returning the correctly renormalized post-measurement statevector over the full register

**Tests** (`tests/`)
- Basic sanity check for `custom_state` + `partial_measurement` (`test_gates.py`)

### Not yet implemented

- **`Circuit` class** (`circuit.py`) — a builder API for composing gates into a named circuit and running it end-to-end, instead of manually threading the statevector through `apply_*` calls
- **Debug utilities** (`debug.py`) — `show_statevector` (pretty-print amplitudes/probabilities) and `check_normalized` (standalone normalization check)
- **`reset_qubit`** — force a qubit back to |0> after measurement
- Broader test coverage — gates, measurement, and edge cases beyond the single existing script
- Packaging cleanup — `LICENSE` file is a placeholder/typo'd as `LISENCE`, `pyproject.toml` author fields still say "Your Name"/"you@example.com", and the GitHub URLs are placeholders

## Design notes

- Gate application comes in two flavors: `apply_gate` builds the full 2^n x 2^n operator via `np.kron` (simple, general, but O(4^n) memory) — used only for single-qubit gates. Multi-qubit gates (`apply_cnot`, `apply_cz`, `apply_swap`, `apply_toffoli`) instead operate directly on state-vector indices using bitmasks, avoiding full matrix construction.
- All qubit indexing is big-endian (qubit 0 is the most significant bit in the basis string).
- `sampling` draws from the probability distribution without collapsing the state; `full_measurement` and `partial_measurement` return properly collapsed and renormalized statevectors.

## License

See `LICENSE` (currently a placeholder — needs to be filled in).
