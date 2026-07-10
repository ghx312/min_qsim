from min_qsim import state, gates, measurement
import numpy as np

vector = [1/np.sqrt(2), 0, 0, 0, 0, 0, 0, 1/np.sqrt(2)]
q = state.custom_state(vector)
q_partial = measurement.partial_measurement(q, 3, [1])

print(q_partial)