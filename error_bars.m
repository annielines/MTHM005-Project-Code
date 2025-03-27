%% plot numerical slices onto contour plot as error bars
clear
syms m rho k_avg p q

% define n+ and n-
n_plus = (1 + m) / 2;
n_minus = (1 - m) / 2;

% define coupled differential equations
dm_dt = 2.*(1-p).*( -n_plus.*(rho./(2.*n_plus)).^q + ...
    n_minus.*(rho./(2.*n_minus)).^q );

drho_dt = (2./k_avg).*( -p.*( n_plus.*(rho./(2.*n_plus)).^q + n_minus.*(rho./(2.*n_minus)).^q ) + ...
    (1-p).*( n_plus .* ( rho./(2*n_plus) ).^q ).*(k_avg - 2.*q - 2.*(k_avg - q).*(rho./(2.*n_plus)) ) + ...
    (1-p).*( n_minus .* ( rho./(2*n_minus) ).^q ).*(k_avg - 2.*q - 2.*(k_avg - q).*(rho./(2.*n_minus)) ) );

% compute jacobian matrix
J = jacobian([dm_dt; drho_dt], [m, rho]);

% define rho*
rho_star = ((1-p)*(k_avg - 2*q) - p) / (2*(1-p)*(k_avg - q));

% evalulate jacobian at steady state (0, rho*)
J_steady_state = subs(J, [m, rho], [0, rho_star]);

% define parameter ranges and average degree
p_values = linspace(0, 0.9, 200);
q_values = linspace(0, 2, 200);
k_value = 6;

% loop over p_values
for i = 1:length(p_values)

    % update current p value
    p_value = p_values(i);

    % loop over q_values
    for j = 1:length(q_values)

        % update current q value
        q_value = q_values(j);

        % evaluate pc
        pc = (k_value - 2*q_value) / (1 + k_value - 2*q_value);

        % if p > pc set to NaN (do not want to plot these)
        if p_value > pc

            max_real_eig_values(i, j) = NaN;

        else

            % substitute p and q to numerically evaluate jacobian
            J_steady_state_num = double(subs(J_steady_state, [p, q, k_avg], [p_value, q_value, k_value]));

            % numerically evaluate eigenvalues of jacobian
            eig_values = eig(J_steady_state_num);

            % take real parts of eigenvalue pair and store maximimum in a matrix
            max_real_eig_values(i, j) = max(real(eig_values));

        end

    end

end

% define contour levels centred around 0
contour_levels = unique([linspace(-2, 0, 8), linspace(0, 0.2, 8)]);
% made these more precise above 0 since (magnitudes of positive max eigenvalues are smaller than negative ones)
% 'unique' so zero is only included once

figure;
hold on
contourf(p_values, q_values, max_real_eig_values', contour_levels, 'LineColor', 'none')
contour(p_values, q_values, max_real_eig_values', contour_levels, 'k', 'LineWidth', 2)
contour(p_values, q_values, max_real_eig_values', [0, 0], 'k', 'LineWidth', 2) % contour line at zero

colormap([linspace(0.05, 0.83, 112)', linspace(0.31, 0.94, 112)', linspace(0.45, 1.00, 112)'; ... % Deeper blue gradient (dark to light)
          linspace(1.00, 0.76, 16)', linspace(0.99, 0.26, 16)', linspace(0.99, 0.25, 16)']); % Red gradient (light to dark)

colorbar;
clim([-2, 0.2])
xlabel('p', 'FontSize', 14);
ylabel('q', 'FontSize', 14);
set(gca, 'FontSize', 16);

% plot p = p_c line
plot1 = plot((k_value - 2.*q_values) ./ (1 + k_value - 2.*q_values), q_values, 'k-.', 'LineWidth', 2);

% the following vector were manually made using 'error_bar_data.mat':

% list of fragmentation points
p_frag = [0.75, 0.75, 0.75, 0.7, 0.7, 0.7, 0.65, 0.65, 0.6, 0.6, 0.6, ...
          0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55];

% list of lower bounds
p_low  = [0.05, 0.05, 0.1, 0.05, 0.05, 0.05, 0.05, 0.1, 0.05, 0.1, 0.1, ...
          0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];

% list of upper bonds
p_high = [0, 0, 0, 0.05, 0, 0, 0.05, 0, 0.05, 0.05, 0, 0.05, 0.05, ...
          0.05, 0.05, 0.05, 0.05, 0.05, 0.1];

% list of q values
q = 0.1:0.1:1.9;

% add error bars to contour plot
plot2 = plot(p_frag, q, 'o', 'MarkerFaceColor', [1, 0.176, 0.161], 'MarkerSize', 8);
hold on

% set height for vertical lines
capHeight = 0.02;

for i = 1:length(p_frag)
    x_left  = p_frag(i) - p_low(i);
    x_right = p_frag(i) + p_high(i);
    q_i = q(i);
    
    % add horizontal error bar
    plot([x_left, x_right], [q_i, q_i], '-', 'color', [1, 0.176, 0.161], 'LineWidth', 2)
    
    % add vertical lines at each end
    plot([x_left, x_left], [q_i - capHeight, q_i + capHeight], '-', 'color', [1, 0.176, 0.161], 'LineWidth', 2.6)
    plot([x_right, x_right], [q_i - capHeight, q_i + capHeight], '-', 'color', [1, 0.176, 0.161], 'LineWidth', 2.6)
end

legend([plot1, plot2], 'Approximated', 'Numerical')
