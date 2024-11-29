clearvars;

L = 600;
alpha_C = 0.168161435;
alpha_F = 0.4506726457;
alpha_A = 0.3811659193;
Beta_F = 0.00005425;
Beta_A = 0.00005425;
M_LA = 308800000;
M_LF = 377000000;
phi_F = 0.00003461169813;
phi_A = 0.00005130434783;
Nu_F = 0.0000000005258143724;
Nu_A = 0.0000000005258143724;
I_F = 345873085 / 100;
I_A = 283377397 / 100;
tau = 0.000000001112433819;

S = 435;
C = 0;
Y = 0;
N = 0;
L_F = 0;
L_A = 0;
M_F = 345873085;
M_A = 283377397;

equation_names = {'S(t)', 'C(t)', 'Y(t)', 'N(t)', 'L_F(t)', 'L_A(t)', 'M_F(t)', 'M_A(t)'};

tspan = [0 100];
dt = 0.001; 

num_steps = ceil((tspan(2) - tspan(1)) / dt);

t = zeros(num_steps + 1, 1);
y = zeros(num_steps + 1, 8);

X = zeros(1000, 8);
Z = zeros(num_steps + 1, 3, 1000);

t(1) = tspan(1);
y(1, :) = [S, C, Y, N, L_F, L_A, M_F, M_A];

for j = 1:1000
    L0 = L - (L / 10) + (2 * L / 10) * rand;
    alpha_C0 = alpha_C - (alpha_C / 10) + (2 * alpha_C / 10) * rand;
    alpha_F0 = alpha_F - (alpha_F / 10) + (2 * alpha_F / 10) * rand;
    alpha_A0 = alpha_A - (alpha_A / 10) + (2 * alpha_A / 10) * rand;
    Beta_F0 = Beta_F - (Beta_F / 10) + (2 * Beta_F / 10) * rand;
    Beta_A0 = Beta_A - (Beta_A / 10) + (2 * Beta_A / 10) * rand;
    M_LA0 = M_LA - (M_LA / 10) + (2 * M_LA / 10) * rand;
    M_LF0 = M_LF - (M_LF / 10) + (2 * M_LF / 10) * rand;
    phi_F0 = phi_F - (phi_F / 10) + (2 * phi_F / 10) * rand;
    phi_A0 = phi_A - (phi_A / 10) + (2 * phi_A / 10) * rand;
    Nu_F0 = Nu_F - (Nu_F / 10) + (2 * Nu_F / 10) * rand;
    Nu_A0 = Nu_A - (Nu_A / 10) + (2 * Nu_A / 10) * rand;
    I_F0 = I_F - (I_F / 10) + (2 * I_F / 10) * rand;
    I_A0 = I_A - (I_A / 10) + (2 * I_A / 10) * rand;
    tau0 = tau - (tau / 10) + (2 * tau / 10) * rand;
    
    for i = 1:num_steps
        current_state = y(i, :);
        dL_F = Nu_F0 * M_LF0 * (L0 - current_state(5) - current_state(6));
        dL_A = Nu_A0 * M_LA0 * (L0 - current_state(6) - current_state(5));
        dM_F = -I_F0;
        dM_A = -I_A0;
        dS = -alpha_C0 * current_state(1) - alpha_F0 * current_state(1) - alpha_A0 * current_state(1);
        dC = alpha_C0 * current_state(1) - Beta_F0 * current_state(2) * current_state(3) - Beta_A0 * current_state(2) * current_state(4) - phi_F0 * current_state(5) * current_state(2) - phi_A0 * current_state(6) * current_state(2) - tau0 * I_F0 * current_state(2) - tau0 * I_A0 * current_state(2);
        dY = alpha_F0 * current_state(1) + Beta_F0 * current_state(2) * current_state(3) + phi_F0 * current_state(5) * current_state(2) + tau0 * I_F0 * current_state(2);
        dN = alpha_A0 * current_state(1) + Beta_A0 * current_state(2) * current_state(4) + phi_A0 * current_state(6) * current_state(2) + tau0 * I_A0 * current_state(2);
    
        y(i + 1, :) = current_state + dt * [dS, dC, dY, dN, dL_F, dL_A, dM_F, dM_A]; 
    
        t(i + 1) = t(i) + dt;
    end
    X(j, :) = y(end, :);
    Z(:, :, j) = y(:, 2:4);
end

for k = 1:8
    subplot(2, 4, k);
    histogram(X(:, k), 20);
    title(['Histogram for ', equation_names{k}]);
    xlabel([equation_names{k}, ' final value']);
    ylabel(['Number of Occurences'])
end

num_lines = 500;

colors = jet(num_lines);
equation_names2 = {'C(t)', 'Y(t)', 'N(t)'};

figure;
for w = 1:3
    subplot(1, 3, w);
    hold on;
    for j = 1:num_lines
        plot(t, Z(:, w, j), 'Color', colors(j, :));
    end
    if w == 1
        plot(t(end), 4, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    end
    if w == 2
        plot(t(end), 219, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    end
    if w == 3
        plot(t(end), 212, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    end
    
    hold off;
    title(['Trajectories for ', equation_names2{w}]);
    xlabel('Time (days)')
    ylabel(equation_names2{w})
end
