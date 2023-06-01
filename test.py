import numpy as np
from build.models import LinearElastic, MCC, SMCC    # Note that this import statement is importing the current build directly from the build folder!
from matplotlib import pyplot as plt

# # LinearElastic test.
# increments = 100
# parameters = np.array([5000.0, 2000.0])
# state = np.array([])
# stress = np.array([-10, -10, -10, 0, 0, 0])
# ea_max = 0.5
# ea_increment = ea_max/increments
# Delta_epsilon_tilde = np.array([-ea_increment, ea_increment/2, ea_increment/2, 0.0, 0.0, 0.0])
# axial_strain = np.arange(0, increments*ea_increment, ea_increment)

# model = LinearElastic(parameters=parameters, state=state, log_severity="verbose")
# print("Model name: {}".format(model.name))
# print("Model type: {}".format(model.type))
# model.set_sigma_prime_tilde(stress)
# p = np.zeros(increments)
# q = np.zeros(increments)
# sigma_prime = np.zeros([increments, 6])
# p[0] = model.p_prime
# q[0] = model.q
# sigma_prime[0,:] = model.sigma_prime_tilde
# print("Increment {}: p = {:.2f} ; q = {:.2f} ; q/p = {:.2f}".format(0, p[0], q[0], q[0]/p[0]))
# i = 0
# while i<increments-1:
#     i += 1
#     model.set_Delta_epsilon_tilde(Delta_epsilon_tildt
#     p[i] = model.p_prime
#     q[i] = model.q
#     sigma_prime[i,:] = model.sigma_prime_tilde
#     print("Increment {}: p = {:.2f} ; q = {:.2f} ; q/p = {:.2f}".format(i, p[i], q[i], q[i]/p[i]))

# plt.plot(axial_strain, q/p)
# plt.xlabel(r"$\epsilon_{a}$ (-)")
# plt.ylabel("q/p (-)")
# plt.show()

# MCC test.
increments = 100
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
state = np.array([1.7477796692480023, 10])
stress = np.array([-10, -10, -10, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([-ea_increment, ea_increment/2, ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = MCC(log_severity="verbose", parameters=parameters, state=state)
print("Model name: {}".format(model.name))
print("Model type: {}".format(model.type))
model.set_sigma_prime_tilde(stress)
p = np.zeros(increments)
q = np.zeros(increments)
sigma_prime = np.zeros([increments, 6])
p[0] = model.p_prime
q[0] = model.q
sigma_prime[0,:] = model.sigma_prime_tilde
print("Increment {}: p = {:.2f} ; q = {:.2f} ; q/p = {:.2f}".format(0, p[0], q[0], q[0]/p[0]))
i = 0
while i<increments-1:
    i += 1
    model.set_Delta_epsilon_tilde(Delta_epsilon_tilde)
    model.solve()
    p[i] = model.p_prime
    q[i] = model.q
    sigma_prime[i,:] = model.sigma_prime_tilde
    print("Increment {}: p = {:.2f} ; q = {:.2f} ; q/p = {:.2f}".format(i, p[i], q[i], q[i]/p[i]))

plt.plot(axial_strain, q)
plt.xlabel(r"$\epsilon_{a}$ (-)")
plt.ylabel("q (kPa)")
plt.show()

# SMCC test.
increments = 100
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02, 4, 0.2, 0.5])
state = np.array([1.7477796692480023, 10, 4])
stress = np.array([-40, -40, -40, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([-ea_increment, ea_increment/2, ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = SMCC(log_severity="verbose", parameters=parameters, state=state)
print("Model name: {}".format(model.name))
print("Model type: {}".format(model.type))
model.set_sigma_prime_tilde(stress)
p = np.zeros(increments)
q = np.zeros(increments)
sigma_prime = np.zeros([increments, 6])
p[0] = model.p_prime
q[0] = model.q
sigma_prime[0,:] = model.sigma_prime_tilde
print("Increment {}: p = {:.2f} ; q = {:.2f} ; q/p = {:.2f}".format(0, p[0], q[0], q[0]/p[0]))
i = 0
while i<increments-1:
    i += 1
    model.set_Delta_epsilon_tilde(Delta_epsilon_tilde)
    model.solve()
    p[i] = model.p_prime
    q[i] = model.q
    sigma_prime[i,:] = model.sigma_prime_tilde
    print("Increment {}: p = {:.2f} ; q = {:.2f} ; q/p = {:.2f}".format(i, p[i], q[i], q[i]/p[i]))
    
plt.plot(axial_strain, q)
plt.xlabel(r"$\epsilon_{a}$ (-)")
plt.ylabel("q (kPa)")
plt.show()