import numpy as np
from geomat.abstract import Elastoplastic
from geomat.utilities import Derivatives
from geomat.models import LinearElastic, MCC, SMCC, C2MC, EMC
from matplotlib import pyplot as plt

class pyMCC(Elastoplastic):

    def __init__(self, parameters, state, log_severity):
        super().__init__()
        self.set_model_name("MCC")
        self.set_model_type("Elastoplastic")

        # Define individual parameters and state variables.
        self.parameters = parameters
        self.state = state
        self.parameters_required = 5
        self.state_required = 1
        self.M = self.parameters[0]
        self.nu = self.parameters[1]
        self.N = self.parameters[2]
        self.lambda_star = self.parameters[3]
        self.kappa_star = self.parameters[4]
        self.p_c = self.state[0]

        # Set the log severity.
        self.initialise_log(log_severity)

        # Check the number of state variables required.
        self.check_inputs(self.name, self.parameters.shape[0], self.state.shape[0], self.parameters_required, self.state_required)

    def compute_D_e(self, sigma_prime, Delta_epsilon):
        # Pressure dependent isotropic elasticity.
        Delta_epsilon_e_vol = self.compute_Delta_epsilon_vol(Delta_epsilon)
        p_prime = self.compute_p_prime(sigma_prime)
        K = self.compute_K_Butterfield(p_prime, Delta_epsilon_e_vol, self.kappa_star, self.settings.EPS)
        G = self.compute_G_given_K_nu(K, self.nu)
        D_e = self.compute_isotropic_linear_elastic_matrix(K, G)
        return D_e

    def compute_f(self, sigma_prime, state):
        # Stress invariants.
        q = self.compute_q(sigma_prime)
        p_prime = self.compute_p_prime(sigma_prime)
        
        # State variables.
        p_c = state[0]
        
        # Yield surface function.
        f = q**2 + self.M**2*p_prime*(p_prime-p_c)

        return f

    def compute_derivatives(self, sigma_prime, state):
        # State variables.
        p_c = state[0]

        # Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
        df_dsigma_prime = np.array((3,3))
        dg_dsigma_prime = np.array((3,3))
        H_s = np.array(np.shape(state))
        B_s = np.array(np.shape(state))
        q = self.compute_q(sigma_prime)
        p_prime = self.compute_p_prime(sigma_prime)
        s = self.compute_s(sigma_prime, p_prime)
        dq_dsigma_prime = self.compute_dq_dsigma_prime(sigma_prime, s, q)
        
        # Yield surface derivatives.
        df_dq = 2*q
        df_dp_prime = self.M**2*(2*p_prime-p_c)
        dp_prime_dsigma_prime = 1.0/3.0*np.eye(3)

        # Derivatives of yield surface and plastic potential function.
        df_dsigma_prime = (df_dp_prime*dp_prime_dsigma_prime) + (df_dq*dq_dsigma_prime)
        dg_dsigma_prime = df_dsigma_prime; # Associated flow.

        # Hardening moduli.
        H_s[0] = (self.M**2*p_prime*p_c)/(self.lambda_star-self.kappa_star)*np.trace(dg_dsigma_prime)

        # State variable increment factors.
        dg_dp_c = -self.M**2*p_prime
        B_s[0] = H_s[0]/dg_dp_c

        # Return Derivatives object.
        derivatives = Derivatives()
        derivatives.df_dsigma_prime = df_dsigma_prime
        derivatives.dg_dsigma_prime = dg_dsigma_prime
        derivatives.H_s = H_s
        derivatives.B_s = B_s
        return derivatives

    def get_state_variables(self):
        return self.state
    
    def set_state_variables(self, state):
        self.state = state

# LinearElastic test.
increments = 100
parameters = np.array([500.0, 1000.0])
state = np.array([])
stress = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
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

# pyMCC test.
increments = 100
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
state = np.array([50.0])
stress = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = pyMCC(parameters=parameters, state=state, log_severity="none")
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

plt.plot(axial_strain, q, "--")
plt.xlabel(r"$\epsilon_{a}$ (-)")
plt.ylabel("q (kPa)")   

# MCC test.
increments = 100
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
state = np.array([50.0])
stress = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = MCC(parameters=parameters, state=state, log_severity="none")
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

plt.plot(axial_strain, q, "-.")
plt.xlabel(r"$\epsilon_{a}$ (-)")
plt.ylabel("q (kPa)")

# SMCC test.
increments = 100
parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02, 5, 1.0, 0.5])
state = np.array([10.0, 5.0])
stress = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
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

# C2MC test.
increments = 100
parameters = np.array([1000.0, 0.0, 0.0, 15.0, 3.0, 29.0, 1.1])
state = np.array([])
stress = np.array([50, 50, 50, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = C2MC(log_severity="none", parameters=parameters, state=state)
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

# EMC test.
increments = 100
parameters = np.array([1000.0, 0.2, 0.0, 40.0, 30.0, 0.0015, 29.0, 1.1])
state = np.array([10.0])
stress = np.array([50, 50, 50, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = EMC(log_severity="none", parameters=parameters, state=state)
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