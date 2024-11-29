import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from pyDOE import lhs
import pandas as pd
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler

def system(y, t, params):
    Beta_F, Beta_A, phi_F, phi_A, tau_F, tau_A = params
    S, C, Y, N, L_F, L_A, M_F, M_A = y
    L = 600
    alpha_C = 0.168161435
    alpha_F = 0.4506726457
    alpha_A = 0.3811659193
    M_LA = 308800000
    M_LF = 377000000
    Nu_F = 0.0000000005258143724
    Nu_A = 0.0000000005258143724
    I_F = 345873.085
    I_A = 283377.397

    dL_F = Nu_F * M_LF * (L - L_F - L_A)
    dL_A = Nu_A * M_LA * (L - L_A - L_F)
    dM_F = -I_F
    dM_A = -I_A
    dS = -alpha_C * S - alpha_F * S - alpha_A * S
    dC = alpha_C * S - Beta_F * C * Y - Beta_A * C * N - phi_F * L_F * C - phi_A * L_A * C - tau_F * I_F * C - tau_A * I_A * C
    dY = alpha_F * S + Beta_F * C * Y + phi_F * L_F * C + tau_F * I_F * C
    dN = alpha_A * S + Beta_A * C * N + phi_A * L_A * C + tau_A * I_A * C

    return [dS, dC, dY, dN, dL_F, dL_A, dM_F, dM_A]

# Function to integrate system and categorize result
def integrate_and_categorize(params, initial_conditions, t):
    sol = odeint(system, initial_conditions, t, args=(params,))
    compartment_Y = sol[:, 2]
    compartment_N = sol[:, 3]
    if compartment_Y[-1] > compartment_N[-1]:
        return (*params, 1)
    elif compartment_Y[-1] < compartment_N[-1]:
        return (*params, 0)
    else:
        return (*params, 2)

# Initial conditions
S0 = 435
M_F0 = 345873085
M_A0 = 283377397
initial_conditions = [S0, 0, 0, 0, 0, 0, M_F0, M_A0]

# Time points
t = np.linspace(0, 1500, int(1500 / 0.001) + 1)

# Parameter bounds based on the given constraints
L = 600
param_bounds = [
    (0, 1/S0),  # Beta_F
    (0, 1/S0),  # Beta_A
    (0, 1/(L + S0)),  # phi_F
    (0, 1/(L + S0)),  # phi_A
    (0, 1/(S0 + M_F0 + M_A0)),  # tau_F
    (0, 1/(S0 + M_F0 + M_A0))   # tau_A
]

# Latin Hypercube Sampling
num_samples = 1000
lhs_samples = lhs(len(param_bounds), samples=num_samples)

# Scale samples to the parameter bounds
param_samples = []
for sample in lhs_samples:
    scaled_sample = [low + (high - low) * value for value, (low, high) in zip(sample, param_bounds)]
    param_samples.append(scaled_sample)

# Use joblib to parallelize the integration and categorization
results = Parallel(n_jobs=-1)(delayed(integrate_and_categorize)(params, initial_conditions, t) for params in param_samples)

# Convert results to numpy array for easier manipulation
results = np.array(results)

# Prepare DataFrame for plotting
df = pd.DataFrame(results, columns=['Beta_F', 'Beta_A', 'phi_F', 'phi_A', 'tau_F', 'tau_A', 'category'])
df['category'] = df['category'].astype(int)
# Filter out the uncertain category
df = df[df['category'] != 2]

# Standardize the features
scaler = StandardScaler()
X = scaler.fit_transform(df[['phi_F', 'phi_A']])
y = df['category']

# Train SVM model
svm_model = SVC(kernel='linear', C=1.0)
svm_model.fit(X, y)
'''
# Plot Parallel Coordinates for all parameters
plt.figure(figsize=(15, 7))
parallel_coordinates(df, 'category', color=('#77AC30', '#FF0000', '#FFA500'), alpha=0.5)
plt.title('Parallel Coordinates Plot of Parameters')
plt.xlabel('Parameters')
plt.ylabel('Values')
plt.show()
'''
scaler = StandardScaler()
W = scaler.fit_transform(df[['Beta_F', 'Beta_A']])
Z = df['category']

# Train SVM model
svm_model = SVC(kernel='linear', C=1.0)
svm_model.fit(W, Z)

labels = {0: "Bill Failed", 1: "Bill Passed", 2: "Uncertain"}
plt.figure(figsize=(10, 5))
for category, color in zip(df['category'].unique(), ['#77AC30', '#A2142F', '#FFA500']):
    subset = df[df['category'] == category]
    plt.scatter(subset['Beta_F'], subset['Beta_A'], color=color, label= labels[category], alpha=0.5)
    ww, zz = np.meshgrid(np.linspace(df['Beta_F'].min(), df['Beta_F'].max(), 100),
                         np.linspace(df['Beta_A'].min(), df['Beta_A'].max(), 100))

    grid = np.c_[ww.ravel(), zz.ravel()]
    grid_scaled = scaler.transform(grid)
    decision_boundary = svm_model.decision_function(grid_scaled).reshape(ww.shape)

    plt.contour(ww, zz, decision_boundary, levels=[0], linewidths=2, colors='black')
plt.xlabel(r'$\beta_F$')
plt.ylabel(r'$\beta_A$')
plt.legend(loc = "upper right")
plt.show()

# Standardize the features
scaler = StandardScaler()
X = scaler.fit_transform(df[['phi_F', 'phi_A']])
y = df['category']

# Train SVM model
svm_model = SVC(kernel='linear', C=1.0)
svm_model.fit(X, y)

labels = {0: "Bill Failed", 1: "Bill Passed", 2: "Uncertain"}
plt.figure(figsize=(10, 5))
for category, color in zip(df['category'].unique(), ['#77AC30', '#A2142F', '#FFA500']):
    subset = df[df['category'] == category]
    plt.scatter(subset['phi_F'], subset['phi_A'], color=color, label= labels[category], alpha=0.5)
    xx, yy = np.meshgrid(np.linspace(df['phi_F'].min(), df['phi_F'].max(), 100),
                         np.linspace(df['phi_A'].min(), df['phi_A'].max(), 100))

    grid = np.c_[xx.ravel(), yy.ravel()]
    grid_scaled = scaler.transform(grid)
    decision_boundary = svm_model.decision_function(grid_scaled).reshape(xx.shape)

    plt.contour(xx, yy, decision_boundary, levels=[0], linewidths=2, colors='black')
plt.xlabel(r'$\phi_F$')
plt.ylabel(r'$\phi_A$')
plt.legend(loc = "upper right")
plt.show()


labels = {0: "Bill Failed", 1: "Bill Passed", 2: "Uncertain"}
plt.figure(figsize=(10, 5))
for category, color in zip(df['category'].unique(), ['#77AC30', '#A2142F', '#FFA500']):
    subset = df[df['category'] == category]
    plt.scatter(subset['tau_F'], subset['tau_A'], color=color, label= labels[category], alpha=0.5)
plt.xlabel(r'$\tau_F$')
plt.ylabel(r'$\tau_A$')
plt.legend(loc = "upper right")
plt.show()
