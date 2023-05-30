import numpy as np
from build.models import MCC    # Note that this import statement is importing the current build directly from the build folder!
from matplotlib import pyplot as plt

increments = 1000
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
state = np.array([1.7477796692480023, 10])
stress = np.array([-10, -10, -10, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([-ea_increment, ea_increment/2, ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = MCC(parameters, state)
model.set_sigma_prime_tilde(stress)
p = np.zeros(increments)
q = np.zeros(increments)
p[0] = model.p_prime
q[0] = model.q

i = 0
# print("Increment {}: p = {} ; q = {} ; q/p = {}".format(i, p[i], q[i], q[i]/p[i]))
while i<increments-1:
    i += 1
    model.set_Delta_epsilon_tilde(Delta_epsilon_tilde)
    model.solve()
    p[i] = model.p_prime
    q[i] = model.q
    # print("Increment {}: p = {} ; q = {} ; q/p = {}".format(i, p[i], q[i], q[i]/p[i]))

plt.plot(axial_strain, q/p)
plt.xlabel(r"$\epsilon_{a}$ (-)")
plt.ylabel("q/p (-)")
plt.show()


    