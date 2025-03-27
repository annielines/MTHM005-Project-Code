%% linear stability analysis for dynamically active steady state
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

            % substitue p and q to numerically evaluate jacobian
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
% made these more precise above 0 since (magnitudes of positive max
% eigenvalues are smaller than negative ones)
% 'unique' so zero is only included once

figure;
hold on
contourf(p_values, q_values, max_real_eig_values', contour_levels, 'LineColor', 'none')
contour(p_values, q_values, max_real_eig_values', contour_levels, 'k', 'LineWidth', 2)
contour(p_values, q_values, max_real_eig_values', [0, 0], 'k', 'LineWidth', 2) % contour line at zero

colormap([linspace(0, 0.4, 112)', linspace(0, 0.4, 112)', linspace(0.3, 0.8, 112)'; ...
    linspace(0.8, 0.3, 16)', linspace(0.3, 0, 16)', linspace(0.3, 0, 16)']);
colorbar;
clim([-2, 0.2])
xlabel('p', 'FontSize', 14);
ylabel('q', 'FontSize', 14);
set(gca, 'FontSize', 16);

% plot p = p_c line
plot1 = plot((k_value - 2.*q_values) ./ (1 + k_value - 2.*q_values), q_values, 'k-.', 'LineWidth', 2);
legend([plot1], 'p = p_c')
