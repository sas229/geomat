import numpy as np
from geomat.models import LinearElastic, MCC, SMCC, Elastoplastic
from matplotlib import pyplot as plt

# LinearElastic test.
increments = 100
parameters = np.array([500.0, 1000.0])
state = np.array([])
stress = np.array([50, 50, 50, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = LinearElastic(parameters=parameters, state=state, log_severity="none")
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

plt.plot(axial_strain, q/p)
plt.xlabel(r"$\epsilon_{a}$ (-)")
plt.ylabel("q/p (-)")
# plt.show()

# MCC test.
increments = 100
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
state = np.array([1.7477796692480023, 50])
stress = np.array([50, 50, 50, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = MCC(log_severity="none", parameters=parameters, state=state)
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
# plt.show()

# SMCC test.
increments = 100
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02, 5, 1.0, 0.5])
state = np.array([1.7477796692480023, 10, 5])
stress = np.array([50, 50, 50, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = SMCC(log_severity="none", parameters=parameters, state=state)
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

delta = 0.25
xrange = np.arange(0, 100+delta, delta)
yrange = np.arange(0, 100+delta, delta)
p, q = np.meshgrid(xrange,yrange)

M = 1
p0 = 100
 
f = q**2+M**2*p*(p-p0)

print(np.shape(f))
plt.contour(f, [0])
plt.show()