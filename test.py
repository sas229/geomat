import numpy as np
from build.models import MCC    # Note that this import statement is importing the current build directly from the build folder!
from matplotlib import pyplot as plt

increments = 100
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
sigma_prime = np.zeros([increments, 6])
p[0] = model.p_prime
q[0] = model.q
sigma_prime[0,:] = model.sigma_prime_tilde
print("Increment {}: p = {} ; q = {} ; q/p = {}".format(0, p[0], q[0], q[0]/p[0]))
print("Increment {}: sigma_11 = {} ; sigma_22 = {} ; sigma_33 = {} ; sigma_12 = {} ; sigma_23 = {} ; sigma_13 = {}".format(0, sigma_prime[0,0], sigma_prime[0,1], sigma_prime[0,2], sigma_prime[0,3], sigma_prime[0,4], sigma_prime[0,5]))
i = 0
while i<increments-1:
    i += 1
    model.set_Delta_epsilon_tilde(Delta_epsilon_tilde)
    model.solve()
    p[i] = model.p_prime
    q[i] = model.q
    sigma_prime[i,:] = model.sigma_prime_tilde
    print("Increment {}: p = {} ; q = {} ; q/p = {}".format(i, p[i], q[i], q[i]/p[i]))
    print("Increment {}: sigma_11 = {} ; sigma_22 = {} ; sigma_33 = {} ; sigma_12 = {} ; sigma_23 = {} ; sigma_13 = {}".format(i, sigma_prime[i,0], sigma_prime[i,1], sigma_prime[i,2], sigma_prime[i,3], sigma_prime[i,4], sigma_prime[i,5]))
    
plt.plot(axial_strain, q/p)
plt.xlabel(r"$\epsilon_{a}$ (-)")
plt.ylabel("q/p (-)")
plt.show()