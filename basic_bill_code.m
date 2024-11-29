% Define the parameters
%L = 600;
alpha_C = 0.168;
alpha_F = 0.449;
alpha_A = 0.382;
Beta_F = 0.00005425;
Beta_A = 0.00005425;
%M_LA = 308800000;
%M_LF = 377000000;
%phi_F = 0.00003461169813;
%phi_A = 0.00005130434783;
%Nu_F = 0.0000000005258143724;
%Nu_A = 0.0000000005258143724;
%I_F = 345873.085;
%I_A = 283377.397;
%tau = 0.000000001112433819;

% Initial state values
S = 435;
C = 0;
Y = 0;
N = 0;
%L_F = 0;
%L_A = 0;
%M_F = 345873085;
%M_A = 283377397;

% Time vector
tspan = [0 100];
dt = 0.001;  % Time step size

% Number of time steps
num_steps = ceil((tspan(2) - tspan(1)) / dt);

% Preallocate arrays to store results
t = zeros(num_steps + 1, 1);
y = zeros(num_steps + 1, 4);

% Initial conditions
t(1) = tspan(1);
y(1, :) = [S, C, Y, N];  %, M_F, M_A, L_F, L_A

% Euler method
for i = 1:num_steps
    % Current state
    current_state = y(i, :);
    
    % Compute derivatives
    %dL_F = Nu_F * M_LF * (L - current_state(5) - current_state(6));
    %dL_A = Nu_A * M_LA * (L - current_state(6) - current_state(5));
    %dM_F = -I_F;
    %dM_A = -I_A;
    dS = -alpha_C * current_state(1) - alpha_F * current_state(1) - alpha_A * current_state(1);
    dC = alpha_C * current_state(1) - Beta_F * current_state(2) * current_state(3) - Beta_A * current_state(2) * current_state(4);% - tau * I_F * current_state(2) - tau * I_A * current_state(2); % - phi_F * current_state(5) * current_state(2) - phi_A * current_state(6) * current_state(2)
    dY = alpha_F * current_state(1) + Beta_F * current_state(2) * current_state(3);% + tau * I_F * current_state(2);%  + phi_F * current_state(5) * current_state(2)
    dN = alpha_A * current_state(1) + Beta_A * current_state(2) * current_state(4);% + tau * I_A * current_state(2); % + phi_A * current_state(6) * current_state(2)
    
    % Update state using Euler method
    y(i + 1, :) = current_state + dt * [dS, dC, dY, dN]; %, dM_F, dM_A, dL_F, dL_A
    
    
    % Update time
    t(i + 1) = t(i) + dt;
end

% Create plot
figure;
hold on;

% Plot each series and assign handles for legend
h1 = plot(t, y(:,1), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);  % S(t)
h2 = plot(t, y(:,2), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);  % C(t)
h3 = plot(t, y(:,3), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);  % Y(t)
h4 = plot(t, y(:,4), 'Color', 'k', 'LineWidth', 2);  % N(t) [0.6350 0.0780 0.1840]

% Plot specific points and assign handles for legend
h7 = plot(70, 219, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.5 0]);  % Point at (70, 219) with green face
h8 = plot(70, 212, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');  % Point at (70, 212) with maroon face [0.6350 0.0780 0.1840]
h9 = plot(70, 4, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', '[0.6350 0.0780 0.1840]');  % Point at (70, 212) with maroon face 
% Set axis labels and title
title({'Predicted Votes for the ACESA in the House',' Without Lobbying or Donations'});
xlabel('Time');
ylabel('Number of People');

% Set x-axis limit
xlim([0 70]);

% Legend
legend([h1, h2, h3, h4, h7, h8, h9], {'S(t)', 'C(t)', 'Y(t)', 'N(t)', 'Actual Yes Vote', 'Actual No Votes', 'Not Voting'}, 'Location', 'bestoutside', 'FontSize', 12);

% Grid
grid on;

% Hold off to finish the plot
hold off;
