!pip install pyDOE2
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from sklearn.preprocessing import scale
from scipy.stats import rankdata
from sklearn.linear_model import LinearRegression
import pyDOE2
import matplotlib.pyplot as plt




def model(y, t, params):
    # Unpack parameters
    alpha_F = params['alpha_F']
    alpha_A = params['alpha_A']
    alpha_C = params['alpha_C']
    nu_FM_LF = params['nu_FM_LF']
    nu_AM_LA = params['nu_AM_LA']
    #M_LF = params['M_LF']
    #M_LA = params['M_LA']
    beta_F = params['beta_F']
    beta_A = params['beta_A']
    phi_F = params['phi_F']
    phi_A = params['phi_A']
    tau_FI_F = params['tau_FI_F']
    tau_AI_A = params['tau_AI_A']
    #I_F = params['I_F']  # Constant value for I_F
    #I_A = params['I_A']  # Constant value for I_A

    S, L_F, L_A, C, Y, N = y #M_F, M_A,

    dS = -(alpha_F + alpha_A + alpha_C) * S
    dL_F = nu_FM_LF * (1 - L_F - L_A)  # Assuming L = 1
    dL_A = nu_AM_LA * (1 - L_F - L_A)  # Assuming L = 1
    #dM_F = -I_F
    #dM_A = -I_A
    dC = alpha_C * S - beta_F * C * Y - beta_A * C * N - phi_F * L_F * C - phi_A * L_A * C - (
                tau_FI_F + tau_AI_A) * C
    dY = alpha_F * S + beta_F * C * Y + phi_F * L_F * C + tau_FI_F * C
    dN = alpha_A * S + beta_A * C * N + phi_A * L_A * C + tau_AI_A * C

    # Pack derivatives
    return [dS, dL_F, dL_A, dC, dY, dN] #dM_F, dM_A,


num_samples = 10000
param_names = ['alpha_F', 'alpha_A', 'alpha_C', 'nu_FM_LF', 'nu_AM_LA', 'beta_F', 'beta_A', 'phi_F', 'phi_A', 'tau_FI_F', 'tau_AI_A']

param_ranges = np.array([
    [0, 1],  # alpha_F
    [0, 1],  # alpha_A
    [0, 1],  # alpha_C
    [0, 1],  # nu_FM_LF
    [0, 1],  # nu_AM_LA
    #[0, 0],  # M_LF
    #[0, 0],  # M_LA
    [0, 1],  # beta_F
    [0, 1],  # beta_A
    [0, 1],  # phi_F
    [0, 1],  # phi_A
    [0, 0.1],  # tau_F
    [0, 0.1],  # tau_A
    #[0.9, 1],  # I_F
    #[0.9, 1]  # I_A
])

param_samples = pyDOE2.lhs(len(param_names), samples=num_samples)
param_samples_scaled = np.zeros(param_samples.shape)
for i in range(len(param_names)):
    param_samples_scaled[:, i] = param_samples[:, i] * (param_ranges[i, 1] - param_ranges[i, 0]) + param_ranges[i, 0]

initial_conditions = [1, 0, 0, 0, 0, 0]  # Example initial conditions
time_span = np.linspace(0, 100, 1000)  # Example time span

outcomes = np.zeros(num_samples)
for i in range(num_samples):
    params = dict(zip(param_names, param_samples_scaled[i]))
    solution = odeint(model, initial_conditions, time_span, args=(params,))

    Y_final = solution[-1, 4]
    N_final = solution[-1, 5]

    outcomes[i] = Y_final - N_final

param_ranks = np.apply_along_axis(rankdata, 0, param_samples_scaled)
outcome_ranks = rankdata(outcomes)

PRCC_values = np.zeros(len(param_names))
for i in range(len(param_names)):
    X = param_ranks[:, i]
    Y = outcome_ranks
    Z = np.delete(param_ranks, i, axis=1)

    reg_X = LinearRegression().fit(Z, X)
    X_res = X - reg_X.predict(Z)

    reg_Y = LinearRegression().fit(Z, Y)
    Y_res = Y - reg_Y.predict(Z)

    PRCC_values[i] = np.corrcoef(X_res, Y_res)[0, 1]

for i in range(len(param_names)):
    print(f'{param_names[i]}: PRCC = {PRCC_values[i]:.4f}')

colors = ['#77AC30' if val > 0 else '#A2142F' for val in PRCC_values]  # Green for > 0, Red for <= 0

param_names_greek = [r'$\alpha_F$', r'$\alpha_A$', r'$\alpha_C$', r'$M_{LF}\nu_{F}$', r'$M_{LA}\nu_{A}$', r'$\beta_F$', r'$\beta_A$', r'$\phi_F$', r'$\phi_A$', r'$I_F\tau_{F}$', r'$I_A\tau_{A}$']

plt.bar(param_names_greek, PRCC_values, color=colors)  # Set the bar color based on the condition
plt.xticks(rotation=45, ha='right')
plt.ylabel('PRCC Value')
plt.title('Parameter Impact on the Success of a Bill')
plt.grid(True)
plt.tight_layout()
plt.show()
