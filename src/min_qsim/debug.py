import numpy as np

def format_amplitude(c: complex, decimals: int=3):
    """
    Formats a complex amplitude for display, rounding to `decimals` places
    and dropping the imaginary part if it rounds to zero.

    Args:
        c: complex amplitude to format
        decimals: number of decimal places to round to

    Returns:
        A clean string like "0.707", "0.5+0.5j", or "-0.707"
    """
    re = round(c.real, decimals)
    im = round(c.imag, decimals)
    if re == 0 and im == 0: return "0"
    elif re == 0: return f"{im}j"
    elif im == 0: return f"{re}"
    else:
        if im < 0: return f"{re}-{abs(im)}j"
        elif im > 0: return f"{re}+{im}j"

def state_to_string(state, n: int, tol:float=1e-6):
    """
    Formats a statevector into a nice string

    Args:
        state: n-qubit system
        n: total number of qubits in the system
        tol: tolerance defaulted to 1e-6

    Returns:
        state_string: statevector formatted into a string
    """
    state_list = []
    for i in range(2 ** n):
        if not np.isclose(abs(state[i]), 0, atol=tol):
            state_list.append(f"{format_amplitude(state[i])}|{bin(i)[2:].zfill(n)}>")
    state_string = " + ".join(i for i in state_list)

    return state_string

