import numpy as np
from models import MCC
from matplotlib import pyplot as plt

increments = 1000
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
state = np.array([1.7477796692480023, 10])
stress = np.array([-10, -10, -10, 0, 0, 0])
delta_epsilon_tilde = np.array([-0.0005, 0.00025, 0.00025, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, 1000*0.0005, 0.0005)

model = MCC(parameters, state)
model.set_sigma_prime(stress)
p = np.zeros(increments)
q = np.zeros(increments)
p[0] = model.p_prime
q[0] = model.q

i = 0
print("Increment {}: p = {} ; q = {} ; q/p = {}".format(i, p[i], q[i], q[i]/p[i]))
while i<=increments-2:
    i += 1
    model.set_strain_increment(delta_epsilon_tilde)
    model.solve()
    p[i] = model.p_prime
    q[i] = model.q
    print("Increment {}: p = {} ; q = {} ; q/p = {}".format(i, p[i], q[i], q[i]/p[i]))

plt.plot(axial_strain, q/p)
plt.show()

    