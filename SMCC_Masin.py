##%% Clear interpreter.
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

# Import modules.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager 
import matplotlib as mpl
import sys

# #%% Set font.
# font_manager.findfont('Times New Roman') 
# plt.rc('font', family='serif') 
# plt.rc('font', serif='Times New Roman') 
# mpl.rcParams['mathtext.fontset'] = 'cm'
# mpl.rcParams['mathtext.rm'] = 'serif'

#%% Subroutines.

# Convert Cauchy stress tensor to Voigt equivalent.
def CauchyToVoigt(S_Cauchy):
    S_Voigt = np.zeros(6)
    S_Voigt[0] = S_Cauchy[0,0]
    S_Voigt[1] = S_Cauchy[1,1]
    S_Voigt[2] = S_Cauchy[2,2]
    S_Voigt[3] = S_Cauchy[2,1]
    S_Voigt[4] = S_Cauchy[2,0]
    S_Voigt[5] = S_Cauchy[1,0]
    return S_Voigt;

#%% Subroutine to convert Voigt stress tensor to Cauchy equivalent.
def VoigtToCauchy(S_Voigt):
    S_Cauchy = np.zeros((3,3))
    S_Cauchy[0,0] = S_Voigt[0]
    S_Cauchy[1,1] = S_Voigt[1]
    S_Cauchy[2,2] = S_Voigt[2]
    S_Cauchy[2,1] = S_Voigt[3]
    S_Cauchy[2,0] = S_Voigt[4]
    S_Cauchy[1,0] = S_Voigt[5]
    S_Cauchy[1,2] = S_Voigt[3]
    S_Cauchy[0,2] = S_Voigt[4]
    S_Cauchy[0,1] = S_Voigt[5]
    return S_Cauchy;

#%% Subroutine to calculate stress invariants.
def StressInvariants(sig):
    # Mean and deviatoric stress tensors (in Cauchy form).
    p = 1/3*np.trace(sig) 
    s = sig-(np.identity(3)*p)
    
    # Stress tensor invariants.
    I1 = np.trace(sig)
    I2 = 1/2*((np.trace(sig))**2-np.trace(sig**2))
    I3 = np.linalg.det(sig)
    
    # Deviatoric stress invariants.
    J1 = 0
    J2 = 1/2*np.trace(s**2)
    J3 = 1/3*np.trace(s**3)
    q = np.sqrt(3*J2)
    return q, p, s, I1, I2, I3, J1, J2, J3;

#%% Subroutine to create (6 x 6) elastic matrix from bulk and shear modulus.
def ElasticMatrix(K,G):
    I_Voigt = CauchyToVoigt(np.identity(3)) # Identity matrix (3 x 3) converted to Voigt notation as outer product is useful in constructing D_e.
    D_e = (K-2/3*G)*np.outer(I_Voigt,I_Voigt)+((2*G)*np.identity(6))
    return D_e;

#%% Subroutine to calculate yield function scalar.
def CheckYieldFunction(q,p,M,p_c,s):
    F = q**2 - M**2*p*(p_c*s-p)
    return F; 

#%% Subroutine to find yield surface intersection using the Pegasus method after Sloan et al. (2001).
def YieldSurfaceIntersection(sig, deps, kap, nu, M, e, p_c, alpha_0, alpha_1, FTOL, MAXITS, debug):
    # Calculate incremental volumetric strain and assume fully elastic.
    deps_v_e = np.trace(VoigtToCauchy(deps))
    q, p, s, I1, I2, I3, J1, J2, J3 = StressInvariants(sig) 
    
    # Calculate secant moduli for alpha = 0 & 1.
    if deps_v_e == 0:
        K_0 = p/kap # Assume tangent bulk modulus.
    else:
        K_0 = p/(alpha_0*deps_v_e)*((np.exp((1+e)*(1-np.exp(alpha_0*deps_v_e)))/kap)-1) # Estimated secant bulk modulus.
    G_0 = (3*(1-(2*nu))/(2*(1+nu)))*K_0
    if deps_v_e == 0:
        K_1 = p/kap # Assume tangent bulk modulus.
    else:
        K_1 = p/(alpha_1*deps_v_e)*((np.exp((1+e)*(1-np.exp(alpha_1*deps_v_e)))/kap)-1) # Estimated secant bulk modulus.
    G_1 = (3*(1-(2*nu))/(2*(1+nu)))*K_1
    
    # Assemble elastic matrices for alpha = 0 & 1.
    D_e_0 = ElasticMatrix(K_0,G_0)
    D_e_1 = ElasticMatrix(K_1,G_1)
    
    # Calculate elastic stress increments for alpha = 0 & 1.
    dsig_0 = alpha_0*np.inner(D_e_0,deps)
    dsig_1 = alpha_1*np.inner(D_e_1,deps)
    
    # Recalculate stress invariants for alpha = 0 & 1.
    sig_0 = sig + VoigtToCauchy(dsig_0)
    sig_1 = sig + VoigtToCauchy(dsig_1)
    q_0, p_0, s_0, I1_0, I2_0, I3_0, J1_0, J2_0, J3_0 = StressInvariants(sig_0) 
    q_1, p_1, s_1, I1_1, I2_1, I3_1, J1_1, J2_1, J3_1 = StressInvariants(sig_1) 
    
    # Check yield surface function values for alpha = 0 & 1.
    F_0 = CheckYieldFunction(q_0, p_0, M, p_c, 1)
    F_1 = CheckYieldFunction(q_1, p_1, M, p_c, 1)
    
    count = 0
    if debug ==1:
        print('Pegasus yield surface intersection algorithm...')
    while count <= MAXITS:
        alpha_n = alpha_1-F_1*(alpha_1-alpha_0)/(F_1-F_0)
        
        # Calculate secoant moduli for new alpha.
        if deps_v_e == 0:
            K_n = p/kap # Assume tangent bulk modulus.
        else:
            K_n = p/(alpha_n*deps_v_e)*((np.exp((1+e)*(1-np.exp(alpha_n*deps_v_e)))/kap)-1) # Estimated secant bulk modulus.
        G_n = (3*(1-(2*nu))/(2*(1+nu)))*K_n
        
        # Assemble elastic matrices for new alpha.
        D_e_n = ElasticMatrix(K_n,G_n)
        
        # Calculate elastic stress increments for new alpha.
        dsig_n = alpha_n*np.inner(D_e_n,deps)
        
        # Recalculate stress invariants for new alpha.
        sig_n = sig + VoigtToCauchy(dsig_n)
        q_n, p_n, s_n, I1_n, I2_n, I3_n, J1_n, J2_n, J3_n = StressInvariants(sig_n) 
        
        # Check yield surface function values for new alpha.
        F_n = CheckYieldFunction(q_n, p_n, M, p_c, 1)
        
        if debug == 1:
            print('Iteration:', count, '; F_0 = ', '%.1e' % F_0, '; F_1 = ', '%.1e' % F_1,'; F_n = ', '%.1e' % F_n, '; alpha =', '%.3f' % alpha_n)
                 
        if F_n*F_1 < 0:
            alpha_0 = alpha_1
            F_0 = F_1
        else:
            alpha_0 = alpha_0
            F_0 = (F_0*F_1)/(F_1+F_n)
        alpha_1 = alpha_n
        F_1 = F_n
        count = count+1
        if np.abs(F_n) <= FTOL:
            break
    return alpha_n;

#%% Subroutine to determine the elastic fraction of an elastoplastic unloading increment using the Pegasus method after Sloan et al. (2001).
def ElastoPlasticUnloading(sig, deps, kap, nu, M, e, p_c, alpha_0, alpha_1, FTOL, NSUB, MAXITS, debug):   
    # Calculate incremental volumetric strain and assume fully elastic.
    deps_v_e = np.trace(VoigtToCauchy(deps))
    q, p, s, I1, I2, I3, J1, J2, J3 = StressInvariants(sig) 
    
    # Loop to find alpha values bounding elastoplastic crossing.
    if debug ==1:
            print('Identify bounding values of alpha for yield surface intersection for elatoplastic unloading...')
    count_1 = 0
    while count_1 < MAXITS:
        # Calculate increment in alpha.
        d_alpha = (alpha_1-alpha_0)/NSUB
        alpha_n = alpha_0+d_alpha
        count_2 = 0
        while count_2 < NSUB:
            # Calculate secant moduli for new alpha.
            if deps_v_e == 0:
                K_n = p/kap # Assume tangent bulk modulus.
            elif alpha_n == 0:
                K_n = p/kap # Assume tangent bulk modulus.
            else:
                K_n = p/(alpha_n*deps_v_e)*((np.exp((1+e[i])*(1-np.exp(alpha_n*deps_v_e)))/kap)-1) # Estimated secant bulk modulus.
            G_n = (3*(1-(2*nu))/(2*(1+nu)))*K_n
            
            # Assemble elastic matrices for new alpha.
            D_e_n = ElasticMatrix(K_n,G_n)
            
            # Calculate elastic stress increments for new alpha.
            dsig_n = alpha_n*np.inner(D_e_n,deps)
            
            # Recalculate stress invariants for new alpha.
            sig_n = sig + VoigtToCauchy(dsig_n)
            q_n, p_n, s_n, I1_n, I2_n, I3_n, J1_n, J2_n, J3_n = StressInvariants(sig_n) 
            
            # Check yield surface function values for new alpha.
            F_n = CheckYieldFunction(q_n, p_n, M, p_c, 1)
            
            if F_n > FTOL:
                alpha_1 = alpha_n
                break
            else:
                alpha_0 = alpha_n
            # Increment alpha and count index.
            alpha_n = alpha_n+d_alpha
            count_2 = count_2+1 
        count_1 = count_1+1
    alpha_0 = alpha_n-d_alpha
    alpha_1 = alpha_n
    if debug ==1:
        print('Lower bound alpha =','%.8f' % alpha_0, '; Upper bound alpha =','%.8f' % alpha_1)
    return alpha_0, alpha_1;
        
#%% Subroutine to calculate derivatives for the MCC model.
def CalculateDerivatives(sig, q, p, M, p_c, e, s, k, A):
    # From yield surface function.
    dF_dp_dp_dsig = M**2*(2*p-p_c*s)*1/3*np.array([1, 1, 1, 0, 0, 0])
    dF_dq_dq_dsig = 3*np.array([(sig[0,0]-p), (sig[1,1]-p), (sig[2,2]-p), 2*sig[2,1], 2*sig[2,0], 2*sig[1,0]])
    dF_dsig = dF_dp_dp_dsig + dF_dq_dq_dsig
    dF_dp_c = -(M**2*p*s)
    dF_ds = -(M**2*p*p_c)
    
    # From plastic potential function.
    dG_dsig = dF_dsig # Associated flow.
    dG_dp = M**2*(2*p-p_c*s)    
    
    # Hardening parameter, H.
    trace_dF_dsig = (dF_dsig[0]+dF_dsig[1]+dF_dsig[2])
    dev_dF_dsig = VoigtToCauchy(dF_dsig)-1/3*np.identity(3)*trace_dF_dsig
    H_p_c = ((M**2*p*p_c)/(lam-kap))*s*trace_dF_dsig
    H_s = ((M**2*p*p_c)/(lam-kap))*-k*(s-1)*np.sqrt((1-A)*trace_dF_dsig**2 + A*2/3*np.tensordot(dev_dF_dsig,dev_dF_dsig))
    H = H_p_c + H_s
    B_p_c = -H_p_c/dF_dp_c
    B_s = -H_s/dF_ds
    return dF_dsig, dF_dp_c, dG_dsig, dG_dp, H, B_p_c, B_s;

#%% Debug.
debug = 0

#%% Initial conditions.
n_inc = 100

# Alternative stress paths to be commented in and out.
ea_max = 0.5
deps = np.array([ea_max/n_inc, -ea_max/(2*n_inc), -ea_max/(2*n_inc), 0, 0, 0,]) # Undrained compression.
# ea_max = 0.1
# deps = np.array([ea_max/n_inc, 0, 0, 0, 0, 0,]) # 1D compression.
# ea_max = 0.01
# deps = np.array([ea_max/n_inc, ea_max/n_inc, ea_max/n_inc, 0, 0, 0,]) # Isotropic compression.

#%% MCC parameters.
lam = 0.08
kap = 0.02
nu = 0.2
M = 0.92
N = 1.195
k = 0.2
A = 1.0

#%% Initialise storage arrays.
e = np.zeros(n_inc)
p_c = np.zeros(n_inc)
q = np.zeros(n_inc)
p = np.zeros(n_inc)
F = np.zeros(n_inc)
eps = np.zeros(n_inc)
s = np.zeros(n_inc)

#%% State variables (assuming isotropic initial conditions).
St = 1
s[0] = St**(lam/(lam-kap))
p_c[0] = 10/s[0]
OCR = 1
sig_0 = np.array([p_c[0]*s[0]/OCR, p_c[0]*s[0]/OCR, p_c[0]*s[0]/OCR, 0, 0, 0]) # Wet of critical trial stresses.
sig_0 = VoigtToCauchy(sig_0)
p_0 = p_c[0]/OCR
p_r = 1
e[0] = np.exp(N-kap*np.log(p_0)-(np.log(p_c[0]/p_r)*(lam-kap)))-1;

# Explicit modified Euler algorithm with substepping after Sloan et al. (2001).
for i in range(10):    
    # Check increment type.
    deps_vol = np.trace(VoigtToCauchy(deps))
    q_0, p_0, d_0, I1_0, I2_0, I3_0, J1_0, J2_0, J3_0 = StressInvariants(sig_0) 
           
    # Check yield surface function values for initial stress.
    F_0 = CheckYieldFunction(q_0, p_0, M, p_c[i], s[i])
    
    # Store stresses if first increment.
    if i == 0:
        q[0] = q_0   
        p[0] = p_0 
        F[0] = F_0
        eps[0] = 0
    else:
        q[i] = q_0
        p[i] = p_0
        F[i] = F_0
        eps[i] = eps[i-1]+deps[0]
    
    # Calculate tangent moduli.
    K_0 = p_0/kap
    G_0 = (3*(1-(2*nu))/(2*(1+nu)))*K_0
    
    # Calculate secant moduli.
    if deps_vol == 0:
        K_0_sec = p_0/kap # Assume tangent bulk modulus.
    else:
        K_0_sec = p_0/deps_vol*((np.exp(1-np.exp(deps_vol))/kap)-1) # Estimated secant bulk modulus.
    G_0_sec = (3*(1-(2*nu))/(2*(1+nu)))*K_0_sec
    
    # Assemble elastic matrices.
    D_e = ElasticMatrix(K_0, G_0)
    D_e_sec = ElasticMatrix(K_0_sec, G_0_sec)
    
    # Calculate elastic stress increments assuming tangent and secant moduli.
    dsig_e = np.inner(D_e, deps)
    dsig_e_sec = np.inner(D_e_sec, deps)
    
    # Recalculate stress invariants.
    sig_e_sec = sig_0 + VoigtToCauchy(dsig_e_sec)
    q_e_sec, p_e_sec, d_e_sec, I1_e_sec, I2_e_sec, I3_e_sec, J1_e_sec, J2_e_sec, J3_e_sec = StressInvariants(sig_e_sec) 
    
    # Check yield surface function values assuming fully elastic increment.
    F_e = CheckYieldFunction(q_e_sec, p_e_sec, M, p_c[i], s[i])
    
    # Caclulate elastic increment fraction, alpha.
    FTOL = 1e-8
    LTOL = 1e-6
    STOL = 1e-4
    EPS = 1e-16
    MAXITS_YSI = 20
    MAXITS_EPT = 3
    MAXITS_YSC = 10
    NSUB = 10
    if F_0 > FTOL:
        print('Error: stress state illegal!')
        sys.exit(0)
    elif F_e <= FTOL:
        # Fully elastic increment...
        alpha = 1
    elif F_0 < -FTOL and F_e > FTOL:
        # Elastic-plastic transition increment...
        alpha = YieldSurfaceIntersection(sig_0, deps, kap, nu, M, e[i], p_c[i], 0, 1, FTOL, MAXITS_YSI, debug)
    elif np.abs(F_0) <= FTOL and F_e > FTOL:
        # Check potential for elastoplastic unloading-reloading...  
        # Check loading condition.
        dF_dsig, dF_dp_c, dG_dsig, dG_dp, H, B_p_c, B_s = CalculateDerivatives(sig_0, q_0, p_0, M, p_c[i], e[i], s[i], k, A)
        cos_theta = (np.inner(dF_dsig,dsig_e))/(np.linalg.norm(dF_dsig,2)*np.linalg.norm(dsig_e,2))
        if cos_theta >= -LTOL:
            # Purely plastic increment...
            alpha = 0
        else:
            # Potential for elastoplastic unloading-reloading...
            alpha_0, alpha_1 = ElastoPlasticUnloading(sig_0, deps, kap, nu, M, e[i], p_c[i], 0, 1, FTOL, NSUB, MAXITS_EPT, debug)
            alpha = YieldSurfaceIntersection(sig_0, deps, kap, nu, M, e[i], p_c[i], alpha_0, alpha_1, FTOL, MAXITS_YSI, debug)

    #%% Update stress state based on elastic increment.
    deps_e = deps*alpha
    deps_e_vol = np.trace(VoigtToCauchy(deps_e))
    q_0, p_0, d_0, I1_0, I2_0, I3_0, J1_0, J2_0, J3_0 = StressInvariants(sig_0) 
    
    # If part of increment is elastic, calculate elastic portion of increment.    
    if alpha > 0:
        # Calculate secant moduli.
        if deps_vol == 0:
            K_0_sec = p_0/kap # Assume tangent bulk modulus.
        else:
            K_0_sec = p_0/deps_e_vol*((np.exp(1-np.exp(deps_e_vol))/kap)-1) # Estimated secant bulk modulus.
        G_0_sec = (3*(1-(2*nu))/(2*(1+nu)))*K_0_sec
        
        # Assemble secant elastic matrix.
        D_e_sec = ElasticMatrix(K_0_sec, G_0_sec)
        
        # Calculate elastic stress increment.
        dsig_e_sec = np.inner(D_e_sec, deps_e)
        
        # Update stress state at onset of plastic flow.
        sig_e = sig_0+VoigtToCauchy(dsig_e_sec)
        
        # Update state variables at onset of plastic flow.
        p_c_e = p_c[i]
        e_e = e[i] - (e[i]*deps_e_vol)
        s_e = s[i] # No softening pre-yield.
    
        # Check yield surface function values at end of elastic fraction of increment.
        q_e, p_e, d_e, I1_e, I2_e, I3_e, J1_e, J2_e, J3_e = StressInvariants(sig_e) 
        F_e = CheckYieldFunction(q_e, p_e, M, p_c_e, s_e)
        
        # Calculate remaining fraction of strain increment.
        deps_p = (1-alpha)*deps
    else:
        deps_p = deps
        sig_e = sig_0
        p_c_e = p_c[i]
        e_e = e[i]
        s_e = s[i]
    
    #%% Substep through plastic strains.
    sig_ep = sig_e
    p_c_ep = p_c_e
    e_ep = e_e
    s_ep = s_e
    ITS_YSC = 0
    if alpha < 1:
        # if debug == 1:
        #     print('Stress integration substepping routine with automatic error control...')
        T = 0
        dT = 1
        dT_min = 1e-6
        n = 0
        EPS = 1e-16
        success = 'No'
        count = 0
        while T < 1.0:
            # Calculate strain increment for increment in pseudo-time dT.
            deps_p_dT = deps_p*dT
            deps_vol_dT = np.trace(VoigtToCauchy(deps_p_dT))
            deps_p_p_dT = deps_p_dT[0]+deps_p_dT[1]+deps_p_dT[2]
            deps_p_q_dT = 1/np.sqrt(3)*np.sqrt(2/3*(((deps_p_dT[0]-deps_p_dT[1])**2)+((deps_p_dT[1]-deps_p_dT[2])**2)+((deps_p_dT[2]-deps_p_dT[0])**2)))
            
            #%% Calculate stress invariants for estimate 1.
            q_1, p_1, d_1, I1_1, I2_1, I3_1, J1_1, J2_1, J3_1 = StressInvariants(sig_ep) 
            
            # Calculate tangent moduli for estimate 1.
            K_1 = p_1/kap
            G_1 = (3*(1-(2*nu))/(2*(1+nu)))*K_1
            
            # Assemble elastic stiffness matrix for estimate 1.
            D_e_1 = ElasticMatrix(K_1, G_1)
            
            # Calculate elastic stress increment for estimate 1.
            dsig_e_1 = np.inner(D_e_1, deps_p_dT)
             
            # Calculate derivatives for stress state for estimate 1.
            dF_dsig_1, dF_dp_c_1, dG_dsig_1, dG_dp_1, H_1, B_p_c_1, B_s_1 = CalculateDerivatives(sig_ep, q_1, p_1, M, p_c_ep, e_ep, s_ep, k, A)

            # Assemble elastoplatic stiffness matrix for estimate 1.
            D_ep_1 = D_e_1 - (np.inner(np.outer(np.inner(D_e_1,dF_dsig_1),dF_dsig_1),D_e_1)/(H_1 + np.inner(np.inner(dF_dsig_1, D_e_1), dG_dsig_1)))
            
            # Calculate increment in plastic multiplier lambda for estimate 1.
            dlambda_1 = np.max(np.inner(dF_dsig_1, dsig_e_1)/(H_1 + np.inner(np.inner(dF_dsig_1, D_e_1), dG_dsig_1)),0)
            
            # Calculate stress increment for estimate 1.
            dsig_1 = np.inner(D_ep_1, deps_p_dT)
            
            # Estimate change in hardening parameter for estimate 1.
            dp_c_1 = dlambda_1*B_p_c_1
            p_c_1 = p_c_ep + dp_c_1
            
            # Estimate change in sensitivity parameter for estimate 1.
            ds_1 = dlambda_1*B_s_1
            s_1 = s_ep + ds_1
            
            # Estimate change in voids ratio for estimate 1.
            de_1 = -(1+e_ep)*deps_vol_dT
            e_1 = e_ep + de_1
            
            # Calculate stress for estimate 1.
            sig_1 = sig_ep + VoigtToCauchy(dsig_1)
            
            #%% Calculate stress invariants for estimate 2.
            q_2, p_2, d_2, I1_2, I2_2, I3_2, J1_2, J2_2, J3_2 = StressInvariants(sig_1) 
            
            # Calculate tangent moduli for estimate 2.
            K_2 = p_2/kap
            G_2 = (3*(1-(2*nu))/(2*(1+nu)))*K_2
            
            # Assemble elastic stiffness matrix for estimate 2.
            D_e_2 = ElasticMatrix(K_2, G_2)
            
            # Calculate elastic stress increment for estimate 2.
            dsig_e_2 = np.inner(D_e_2, deps_p_dT)
             
            # Calculate derivatives for stress state for estimate 2.
            dF_dsig_2, dF_dp_c_2, dG_dsig_2, dG_dp_2, H_2, B_p_c_2, B_s_2 = CalculateDerivatives(sig_1, q_2, p_2, M, p_c_1, e_1, s_1, k, A)
     
            # Assemble elastoplatic stiffness matrix for estimate 2.
            D_ep_2 = D_e_2 - np.inner(np.outer(np.inner(D_e_2,dF_dsig_2),dF_dsig_2),D_e_2)/(H_2 + np.inner(np.inner(dF_dsig_2, D_e_2), dG_dsig_2))
                  
            # Calculate increment in plastic multiplier lambda for estimate 2.
            dlambda_2 = np.max(np.inner(dF_dsig_2,dsig_e_2)/(H_2 + np.inner(np.inner(dF_dsig_2,D_e_2),dG_dsig_2)),0)
            
            # Calculate stress increment for estimate 2.
            dsig_2 = np.inner(D_ep_2, deps_p_dT)
            
            # Estimate change in hardening parameter for estimate 2.
            dp_c_2 = dlambda_2*B_p_c_2
            p_c_2 = p_c_ep + dp_c_2
            
            # Estimate change in sensitivity parameter for estimate 2.
            ds_2 = dlambda_2*B_s_2
            s_2 = s_ep + ds_2
            
            # Estimate change in voids ratio for estimate 2.
            de_2 = -(1+e_ep)*deps_vol_dT
            e_2 = e_ep + de_2
            
            #%% Calculate modified Euler stresses.
            sig_ep_ini = sig_ep + 1/2*(VoigtToCauchy(dsig_1+dsig_2))
            p_c_ep_ini = p_c_ep + 1/2*(dp_c_1+dp_c_2)
            e_ep_ini = e_ep + 1/2*(de_1+de_2)
            s_ep_ini = s_ep + 1/2*(ds_1+ds_2)
            
            #%% Calculate error estimate.
            R_n = 1/2*max((np.linalg.norm(VoigtToCauchy(dsig_2-dsig_1))/(np.linalg.norm(sig_ep_ini))),(np.abs(dp_c_2-dp_c_1)/p_c_ep_ini),(np.abs(ds_2-ds_1)/s_ep_ini),(np.abs(de_2-de_1)/e_ep_ini),EPS)
            if R_n < STOL:  
                # Accept modified Euler estimates.
                sig_uncorr = sig_ep_ini
                p_c_uncorr = p_c_ep_ini
                e_uncorr = e_ep_ini
                s_uncorr = s_ep_ini
                
                # Calculate new step size.
                q_step = min(0.9*np.sqrt(STOL/R_n),1.1)
                
                # Correct stresses and state variables back to yield surface.
                while ITS_YSC <= MAXITS_YSC:
                    # Calculate stress invariants for uncorrected stress state.
                    q_uncorr, p_uncorr, d_uncorr, I1_uncorr, I2_uncorr, I3_uncorr, J1_uncorr, J2_uncorr, J3_uncorr = StressInvariants(sig_uncorr) 
        
                    # Calculate tangent moduli for uncorrected stress state and hardening parameter.
                    K_uncorr = p_uncorr/kap
                    G_uncorr = (3*(1-(2*nu))/(2*(1+nu)))*K_uncorr
                
                    # Assemble elastic stiffness matrix for uncorrected stress state and hardening parameter.
                    D_e_uncorr = ElasticMatrix(K_uncorr, G_uncorr)
                    
                    # Check yield function for uncorrected stress state and hardening parameter.
                    F_uncorr = CheckYieldFunction(q_uncorr, p_uncorr, M, p_c_uncorr, s_uncorr)
                    
                    # If initial yield function value is less than tolerance accept the substep.
                    if np.abs(F_uncorr) <= FTOL:
                        sig_ep = sig_uncorr
                        p_c_ep = p_c_uncorr
                        e_ep = e_uncorr
                        s_ep = s_uncorr
                        break
                    elif np.abs(F_uncorr) > FTOL:
                        # Calculate the derivatives for the uncorrected stress state and hardening parameter.
                        dF_dsig_uncorr, dF_dp_c_uncorr, dG_dsig_uncorr, dG_dp_uncorr, H_uncorr, B_p_c_uncorr, B_s_uncorr = CalculateDerivatives(sig_uncorr, q_uncorr, p_uncorr, M, p_c_uncorr, e_uncorr, s_uncorr, k, A)
                        
                        # Calculate correction factor.
                        delta_lambda = F_uncorr/(H_uncorr + np.inner(np.inner(dF_dsig_uncorr, D_e_uncorr), dG_dsig_uncorr))
                        delta_sig = -delta_lambda*np.inner(D_e_uncorr, dG_dsig_uncorr)
                        delta_p_c = delta_lambda*B_p_c_uncorr
                        delta_e = 0
                        delta_s = 0
                        
                        # Correct stresses and hardening parameter.
                        sig_corr = sig_uncorr + VoigtToCauchy(delta_sig)
                        p_c_corr = p_c_uncorr + delta_p_c
                        e_corr = e_uncorr + delta_e
                        s_corr = s_uncorr + delta_s
                        
                        # Calculate stress invariants for corrected stress state.
                        q_corr, p_corr, d_corr, I1_corr, I2_corr, I3_corr, J1_corr, J2_corr, J3_corr = StressInvariants(sig_corr) 
                        
                        # Calculate yield surface function for corrected stress state.
                        F_corr = CheckYieldFunction(q_corr, p_corr, M, p_c_corr, s_corr)
                        
                        # If correction fails use normal correction.
                        if np.abs(F_corr) > np.abs(F_uncorr):
                            delta_lambda = F_uncorr/(np.inner(dF_dsig, dF_dsig))
                            delta_sig = -delta_lambda*dF_dsig
                            delta_p_c = 0
                            delta_e = 0
                            delta_s = 0

                            # Correct stresses and hardening parameter.
                            sig_corr = sig_uncorr + VoigtToCauchy(delta_sig)
                            p_c_corr = p_c_uncorr + delta_p_c
                            e_corr = e_uncorr + delta_e
                            s_corr = s_uncorr + delta_s

                            # Calculate stress invariants for corrected stress state.
                            q_corr, p_corr, d_corr, I1_corr, I2_corr, I3_corr, J1_corr, J2_corr, J3_corr = StressInvariants(sig_corr) 

                            # Calculate yield surface function for corrected stress state.
                            F_corr = CheckYieldFunction(q_corr, p_corr, M, p_c_corr, s_corr)
                            print("Applied normal correction.")
                          
                        # Update stresses and hardening parameter for next increment.
                        sig_uncorr = sig_uncorr + VoigtToCauchy(delta_sig)
                        p_c_uncorr = p_c_uncorr + delta_p_c
                        e_uncorr = e_uncorr + delta_e
                        s_uncorr = s_uncorr + delta_s
                        
                        ITS_YSC += 1
                        print("ITS_YSC = ", ITS_YSC) 
                        print("F_0 = ", F_0) 
                        print("delta_lambda = ", delta_lambda)
                        print("delta_sig = ", delta_sig)
                        print("delta_p_c = ", delta_p_c)
                        print("dF_dsig_uncorr = ", dF_dsig_uncorr)
                        print("dG_dsig_uncorr = ", dG_dsig_uncorr)
                        print("H = ", H)
                        # print("D_e_uncorr = ", D_e_uncorr)
                        print("R_n = ", R_n)
                        print("F_corr = ", F_corr)
                        print("sig_corr = ", sig_corr)
                        input("")
                    
                    # Break from loop when error is less than the specified tolerance.
                    if np.abs(F_corr) <= FTOL:
                        sig_ep = sig_corr
                        p_c_ep = p_c_corr
                        e_ep = e_corr
                        s_ep = s_corr
                        break
                    # Else if maximum number of yield surface correction iterations reached stop.
                    elif ITS_YSC > MAXITS_YSC and   np.abs(F_corr) > FTOL:
                        print('Yield surface drift correction unsuccessful in increment', i+1, 'substep', n, 'after', ITS_YSC, 'iterative corrections!') 
                        print('|F| =', '%.3e' % np.abs(F_corr), '> |FTOL| =', '%.3e' % FTOL)
                        sys.exit(0)
                
                # If previous step failed limit step size growth.
                if success == 'No':
                    q_step = min(q_step,1)
                success = 'Y'
                
                # Update substep time dT and pseudo-time T.
                dT = q_step*dT       
            else:    
                # Calculate new step size.
                q_step = max(0.9*np.sqrt(STOL/R_n),0.1)
                
                # Update dT.
                dT = q_step*dT
            
            #%% Update pseudo-time T and check (i) that the next step size is not smaller than the minimum size specified and (ii) does not exceed T = 1.
            dT = max(dT, dT_min)
            dT = min(dT, 1-T)

            # If current substep was successful update the pseudo-time T.
            if success == 'Y':
                T = T + dT
            n = n + 1 
                    
        if debug == 1:
            print('Increment:', i+1, '; Substeps:', n, '; R_n =','%.2e' % R_n, '; Success:', success, '; Corrections:', ITS_YSC)
        
        # Update state variables and stress state.
        if i < n_inc-1:
            p_c[i+1] = p_c_ep
            s[i+1] = s_ep
            e[i+1] = e_ep
    else:
        # Update state variables and stress state.
        if i < n_inc-1:
            p_c[i+1] = p_c[i] # No change in hardening parameter as fully elastic increment.
            s[i+1] = s[i] # No change in softening parameter as fully elastic increment.
            e[i+1] = e[i] - (1+e[i])*deps_vol
        
    #%% Update stress state for next increment.
    sig_0 = sig_ep
    
#%% Plot output.
plt.close("all")
f1, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2 , 2, sharex=False)
figure = plt.gcf()
# figure.set_size_inches(8, 6)

ax1.plot(p[0:i],q[0:i], 'b')
ax1.plot([0, 1.2*p_c[0]*s[0]], [0, 1.2*M*p_c[0]*s[0]], '-k')
ax1.plot(np.linspace(0,p_c[0]*s[0],1000), np.sqrt(M**2*np.linspace(0, p_c[0]*s[0], 1000)*(p_c[0]*s[0]-np.linspace(0, p_c[0]*s[0], 1000))), 'r', linewidth=0.5)
ax1.plot(np.linspace(0,p_c[-1]*s[-1],1000), np.sqrt(M**2*np.linspace(0, p_c[-1]*s[-1], 1000)*(p_c[-1]*s[-1]-np.linspace(0, p_c[-1]*s[-1], 1000))), 'g', linewidth=0.5)
ax1.set_xlabel(r'$p^\prime$ (kPa)')
ax1.set_ylabel(r'$q$ (kPa)')
ax1.set_xlim(0,int(max(1.2*p_c[0]*s[0],1.2*p_c[-1]*s[-1])))
ax1.set_ylim(0,int(max(1.2*p_c[0]*s[0],1.2*p_c[-1]*s[-1])))

ax2.semilogx(p[0:i],np.log(1+e[0:i]), 'b')
ax2.semilogx(np.linspace(0.1,1000,1000), N-lam*np.log(np.linspace(0.1,1000,1000)), 'k')
ax2.semilogx(np.linspace(0.1,1000,1000), N+lam*np.log(s[0])-lam*np.log(np.linspace(0.1,1000,1000)), 'r', linewidth=0.5)
ax2.semilogx(np.linspace(0.1,1000,1000), N+lam*np.log(s[-1])-lam*np.log(np.linspace(0.1,1000,1000)), 'g', linewidth=0.5)
ax2.set_xlabel(r'$p^\prime$ (kPa)')
ax2.set_ylabel(r'$ln(1+e)$ (-)')
ax2.set_xlim(1,1000)
ax2.set_ylim(0.5,1.2)

ax3.plot(eps[0:i],s[0:i], 'b')
ax3.set_xlabel(r'$\epsilon_{a}$ (-)')
ax3.set_ylabel(r'$s$ (-)')
ax3.set_xlim(0,ea_max)
ax3.set_ylim(0,1.2*s[0])

ax4.plot(eps[0:i],p_c[0:i], 'b')
ax4.set_xlabel(r'$\epsilon_{a}$ (-)')
ax4.set_ylabel(r'$p_{c}$ (kPa)')
ax4.set_xlim(0,ea_max)
ax4.set_ylim(0,1.2*p_c[-1])

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(top=0.9)
plt.show()

