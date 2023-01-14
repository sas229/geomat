import numpy as np
from models import MCC
# from models import SMCC
from matplotlib import pyplot as plt

increments = 100
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
state = np.array([1.7477796692480023, 10])
# parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02, 3.0, 0.2, 0.5])
# state = np.array([1.7477796692480023, 10, 3.0])
stress = np.array([-9.86, -10, -10, 0, 0, 0])
# stress = np.array([-30, -30, -30, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([-ea_increment, ea_increment/2, ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

# model = MCC(parameters, state)
model = MCC(parameters, state, "verbose")
# model = SMCC(parameters, state)
model.set_sigma_prime_tilde(stress)
model.set_Delta_epsilon_tilde(Delta_epsilon_tilde)
p = np.zeros(increments)
q = np.zeros(increments)
p[0] = model.p_prime
q[0] = model.q
# print("Model name: ", model.name)
# print("Model type: ", model.model_type)
# print("IP number: ", model.IP_number)
# print("Mises stress: ", model.mises_stress)
# print("Maximum shear stress: ", model.max_shear)

i = 0
print("Increment {}: q = {} ; p = {} ; q/p = {}".format(i+1, q[i], p[i], q[i]/p[i]))
while i<10:
    i += 1
    model.solve()
    model.set_Delta_epsilon_tilde(Delta_epsilon_tilde)
    p[i] = model.p_prime
    q[i] = model.q
    print("Increment {}: q = {} ; p = {} ; q/p = {}".format(i+1, q[i], p[i], q[i]/p[i]))

plt.plot(axial_strain, q/p)
plt.xlabel(r"$\epsilon_{a}$ (-)")
plt.ylabel("q/p (-)")
plt.show()


    