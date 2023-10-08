import numpy as np
from geomat.library import LinearElastic, SMCC, C2MC, EMC, Elastoplastic, Derivatives
import geomat.checks as checks
import geomat.logging as logging
from matplotlib import pyplot as plt

class MCC(Elastoplastic):

    def __init__(self, parameters, state, log_severity):
        super().__init__()
        self.parameters = parameters
        self.state = state
        self.log_severity = log_severity
        self.settings = self.get_settings()

        # Define individual parameters and state variables.
        self.M = parameters[0]
        self.nu = parameters[1]
        self.N = parameters[2]
        self.lambda_star = parameters[3]
        self.kappa_star = parameters[4]
        self.p_c = state[0]

        # Set model name and type.
        self.name = "MCC"
        self.type = "Elastoplastic"

        # Set the log severity.
        logging.initialise_log(self.log_severity)

        # Required number of parameters and state variables.
        self.parameters_required = 5
        self.state_required = 1
        
        # Check the number of state variables required.
        checks.check_inputs(self.name, self.parameters.shape[0], self.state.shape[0], self.parameters_required, self.state_required)

    def compute_D_e(self, sigma_prime, Delta_epsilon):
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
        print("Getting state variables...")
        return self.state
    
    def set_state_variables(self, state):
        print("Setting state variables...")
        self.state = state

# parameters = np.array([0.92, 0.2, 1.195, 0.08, 0.02])
# state = np.array([50.0])
# model = MCC(parameters, state, log_severity="verbose")
# sigma_prime_tilde = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
# ea_max = 0.5
# increments = 100
# ea_increment = ea_max/increments
# Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])

# model.set_sigma_prime_tilde(sigma_prime_tilde)
# model.set_Delta_epsilon_tilde(Delta_epsilon_tilde)
# model.solve()
# print("Solved increment.")

# LinearElastic test.
increments = 100
parameters = np.array([500.0, 1000.0])
state = np.array([])
stress = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = LinearElastic(parameters=parameters, state=state, log_severity="info")
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
state = np.array([50.0])
stress = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = MCC(parameters=parameters, state=state, log_severity="verbose")
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
state = np.array([10.0, 5.0])
stress = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
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

# C2MC test.
increments = 100
parameters = np.array([1000.0, 0.0, 0.0, 15.0, 3.0, 29.0, 1.1])
state = np.array([])
stress = np.array([50, 50, 50, 0, 0, 0])
ea_max = 0.5
ea_increment = ea_max/increments
Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
axial_strain = np.arange(0, increments*ea_increment, ea_increment)

model = C2MC(log_severity="info", parameters=parameters, state=state)
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

# # EMC test.
# increments = 100
# parameters = np.array([1000.0, 0.2, 0.0, 40.0, 30.0, 0.0015, 29.0, 1.1])
# state = np.array([10.0])
# stress = np.array([50, 50, 50, 0, 0, 0])
# ea_max = 0.5
# ea_increment = ea_max/increments
# Delta_epsilon_tilde = np.array([ea_increment, -ea_increment/2, -ea_increment/2, 0.0, 0.0, 0.0])
# axial_strain = np.arange(0, increments*ea_increment, ea_increment)

# model = EMC(log_severity="info", parameters=parameters, state=state)
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
#     model.set_Delta_epsilon_tilde(Delta_epsilon_tilde)
#     model.solve()
#     p[i] = model.p_prime
#     q[i] = model.q
#     sigma_prime[i,:] = model.sigma_prime_tilde
#     print("Increment {}: p = {:.2f} ; q = {:.2f} ; q/p = {:.2f}".format(i, p[i], q[i], q[i]/p[i]))
    
# plt.plot(p/p[0], q/p)
# plt.xlabel("p/p_i (-)")
# plt.ylabel("q/p (-)")
# plt.show()
