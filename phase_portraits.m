%% phase portaits in (m, rho)-plane
clear

% define coupled differential equations as functions

% magnetism
dm_dt = @(p, q, rho, n_plus, n_minus) 2.*(1-p).*( -n_plus.*(rho./(2.*n_plus)).^q + ...
    n_minus.*(rho./(2.*n_minus)).^q );

% active edge density
drho_dt = @(p, q, rho, n_plus, n_minus, k_avg) (2./k_avg).*( -p.*( n_plus.*(rho./(2.*n_plus)).^q + n_minus.*(rho./(2.*n_minus)).^q ) + ...
    (1-p).*( n_plus .* ( rho./(2*n_plus) ).^q ).*(k_avg - 2.*q - 2.*(k_avg - q).*(rho./(2.*n_plus)) ) + ...
    (1-p).*( n_minus .* ( rho./(2*n_minus) ).^q ).*(k_avg - 2.*q - 2.*(k_avg - q).*(rho./(2.*n_minus)) ) );

% initialise variables
p_vec = [0.2, 0.78, 0.9];
q_vec = [0.5, 1, 2];
k_avg = 8;

% initialise counter for subplot placement
counter = 0;

% loop over q values
for i = 1:length(q_vec)

    % update current q value
    q = q_vec(i);

    % loop over p values
    for j = 1:length(p_vec)

        % update current p value
        p = p_vec(j);

        % define mesh of m and rho values
        [m, rho] = meshgrid(-1:0.05:1, 0.001:0.05:0.501);

        % evaluate dm/dt and drho/dt vectors over mesh for current p and q
        dm = dm_dt(p, q, rho, 0.5.*(1+m), 0.5.*(1-m));
        drho = drho_dt(p, q, rho, 0.5.*(1+m), 0.5.*(1-m), k_avg);

        % calculate magnitudes of vectors
        magnitude = sqrt(dm.^2 + drho.^2);

        % normalize to unit vectors
        dm_unit = dm ./ magnitude; % m-direction
        drho_unit = drho ./ magnitude; % rho-direction

        % set desired length for all arrows
        arrow_length = 4; % can adjust size here
        dm_plot = dm_unit * arrow_length;
        drho_plot = drho_unit * arrow_length;

        % evaluate rho* for current p and q
        rho_star = ( (1 - p)*(k_avg - 2*q) - p ) / ( 2*(1-p)*(k_avg - q) );

        % update counter
        counter = counter + 1;

        % plot the arrows
        figure(1);
        set(gca,'fontsize',16)
        subplot(3,3,counter)
        hold on
        quiver(m, rho, dm_plot, drho_plot, 'LineWidth', 1.4, 'Color','#AF5F5A');
        xlabel('m');
        ylabel('\rho');

        if (p == p_vec(1)) && (q == q_vec(1))

            % plot predicted stabilities of steady states
            xline(0, 'k', 'LineWidth', 2)
            plot(-1,0, 'ko', 'MarkerSize', 16, 'LineWidth', 2)
            plot(1,0, 'ko', 'MarkerSize', 16, 'LineWidth', 2)
            plot(0,rho_star, 'k.', 'MarkerSize', 50)
            plot(0,0, 'ko', 'MarkerSize', 16, 'LineWidth', 2)

            ode_fun = @(t,y) [ dm_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1))); drho_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1)), k_avg)];

            % plot separatrices
            initial_condition = [-0.99, 0.01]; % (m, rho)
            tspan = [0 10];
            [t,y] = ode45(ode_fun,tspan,initial_condition);
            plot(y(:,1),y(:,2), 'k', 'LineWidth', 2)

            initial_condition = [0.99, 0.01]; % (m, rho)
            tspan = [0 10];
            [t,y] = ode45(ode_fun,tspan,initial_condition);
            plot(y(:,1),y(:,2), 'k', 'LineWidth', 2)
            axis([-1 1 0 0.5])
            title(sprintf('p1 < pc, q = %.2f', q))

        elseif (p == p_vec(1)) && (q == q_vec(2))

            % plot predicted stabilities of steady states
            xline(0, 'k', 'LineWidth', 2)

            ode_fun = @(t,y) [ dm_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1))); drho_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1)), k_avg)];

            % run model and plot end points
            initial_conditions = [linspace(-0.999,0.999,100);linspace(0,0.5,100)]; % (m, rho)
            tspan = [0 100];

            for i = 1:length(initial_conditions)
                initial_condition = [initial_conditions(1, i), initial_conditions(2, i)];
                [t,y] = ode45(ode_fun,tspan,initial_condition);
                y_plot(i, 1:2) = [y(end, 1), y(end, 2)];
            end

            plot(y_plot(:,1),y_plot(:,2), 'k', 'LineWidth', 2)
            axis([-1 1 0 0.5])
            title(sprintf('p1 < pc, q = %.2f', q))

        elseif (p == p_vec(1)) && (q == q_vec(3))

            % plot predicted stabilities of steady states
            xline(0, 'k', 'LineWidth', 2)
            plot(-1,0, 'k.', 'MarkerSize', 50)
            plot(1,0, 'k.', 'MarkerSize', 50)
            plot(0,rho_star, 'ko', 'MarkerSize', 16, 'LineWidth', 2)
            plot(0,0, 'ko', 'MarkerSize', 16, 'LineWidth', 2)

            ode_fun = @(t,y) [ dm_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1))); drho_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1)), k_avg)];

            % plot separatrices
            initial_condition = [-0.01, 0.3125]; % (m, rho)
            tspan = [0 100];
            [t,y] = ode45(ode_fun,tspan,initial_condition);
            plot(y(:,1),y(:,2), 'k', 'LineWidth', 2)

            initial_condition = [0.01, 0.3125]; % (m, rho)
            tspan = [0 100];
            [t,y] = ode45(ode_fun,tspan,initial_condition);
            plot(y(:,1),y(:,2), 'k', 'LineWidth', 2)
            axis([-1 1 0 0.5])
            title(sprintf('p1 < pc, q = %.2f', q))

        elseif (p == p_vec(2)) && (q == q_vec(1))

            % plot predicted stabilities of steady states
            xline(0, 'k', 'LineWidth', 2)
            plot(-1,0, 'ko', 'MarkerSize', 16, 'LineWidth', 2)
            plot(1,0, 'ko', 'MarkerSize', 16, 'LineWidth', 2)
            plot(0,rho_star, 'k.', 'MarkerSize', 50)
            plot(0,0, 'ko', 'MarkerSize', 16, 'LineWidth', 2)

            ode_fun = @(t,y) [ dm_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1))); drho_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1)), k_avg)];

            % plot separatrices
            initial_condition = [-0.99, 0.01]; % (m, rho)
            tspan = [0 4000];
            [t,y] = ode45(ode_fun,tspan,initial_condition);
            plot(y(:,1),y(:,2), 'k', 'LineWidth', 2)

            initial_condition = [0.99, 0.01]; % (m, rho)
            tspan = [0 4000];
            [t,y] = ode45(ode_fun,tspan,initial_condition);
            plot(y(:,1),y(:,2), 'k', 'LineWidth', 2)
            axis([-1 1 0 0.5])
            title(sprintf('p2 < pc, q = %.2f', q))

        elseif (p == p_vec(2)) && (q == q_vec(2))

            % plot predicted stabilities of steady states
            xline(0, 'k', 'LineWidth', 2)

            ode_fun = @(t,y) [ dm_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1))); drho_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1)), k_avg)];

            % run model and plot end points
            initial_conditions = [linspace(-0.999,0.999,100);linspace(0,0.5,100)]; % (m, rho)
            tspan = [0 100];

            for i = 1:length(initial_conditions)
                initial_condition = [initial_conditions(1, i), initial_conditions(2, i)];
                [t,y] = ode45(ode_fun,tspan,initial_condition);
                y_plot(i, 1:2) = [y(end, 1), y(end, 2)];
            end

            plot(y_plot(:,1),y_plot(:,2), 'k', 'LineWidth', 2)
            axis([-1 1 0 0.5])
            title(sprintf('p2 < pc, q = %.2f', q))

        elseif (p == p_vec(2)) && (q == q_vec(3))

            % plot predicted stabilities of steady states
            xline(0, 'k', 'LineWidth', 2)
            plot(-1,0, 'k.', 'MarkerSize', 50)
            plot(1,0, 'k.', 'MarkerSize', 50)
            plot(0,rho_star, 'ko', 'MarkerSize', 16, 'LineWidth', 2)
            plot(0,0, 'ko', 'MarkerSize', 16, 'LineWidth', 2)

            ode_fun = @(t,y) [ dm_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1))); drho_dt(p, q, y(2), 0.5.*(1+y(1)), 0.5.*(1-y(1)), k_avg)];

            % plot separatrices
            initial_condition = [-0.01, 0.0379]; % (m, rho)
            tspan = [0 40000];
            [t,y] = ode45(ode_fun,tspan,initial_condition);
            plot(y(:,1),y(:,2), 'k', 'LineWidth', 2)

            initial_condition = [0.01, 0.0379]; % (m, rho)
            tspan = [0 40000];
            [t,y] = ode45(ode_fun,tspan,initial_condition);
            plot(y(:,1),y(:,2), 'k', 'LineWidth', 2)
            axis([-1 1 0 0.5])
            title(sprintf('p2 < pc, q = %.2f', q))

        elseif ( (p == p_vec(3)) && (q == q_vec(1)) ) || ( (p == p_vec(3)) && (q == q_vec(2)) ) || ( (p == p_vec(3)) && (q == q_vec(3)) )

            % plot predicted stabilities of steady states
            xline(0, 'k', 'LineWidth', 2)
            plot(0,0, 'k.', 'MarkerSize', 50)
            axis([-1 1 0 0.5])
            title(sprintf('p3 > pc, q = %.2f', q))

        end
    end
end

set(gca,'fontsize',16)


%% phase portait in (m, rho)-plane for p = pc case
clear

% define coupled differential equations
dm_dt = @(p, q, rho, n_plus, n_minus) 2.*(1-p).*( -n_plus.*(rho./(2.*n_plus)).^q + ...
    n_minus.*(rho./(2.*n_minus)).^q );

drho_dt = @(p, q, rho, n_plus, n_minus, k_avg) (2./k_avg).*( -p.*( n_plus.*(rho./(2.*n_plus)).^q + n_minus.*(rho./(2.*n_minus)).^q ) + ...
    (1-p).*( n_plus .* ( rho./(2*n_plus) ).^q ).*(k_avg - 2.*q - 2.*(k_avg - q).*(rho./(2.*n_plus)) ) + ...
    (1-p).*( n_minus .* ( rho./(2*n_minus) ).^q ).*(k_avg - 2.*q - 2.*(k_avg - q).*(rho./(2.*n_minus)) ) );

% initialise variables
q_vec = [0.5, 1, 2];
k_avg = 8;

% initialise counter for subplot placement
counter = 0;

% loop over q values
for i = 1:length(q_vec)

    % update current q value
    q = q_vec(i);

    % evaluate corresponding p value, where p = pc
    p = (k_avg - 2*q) / (1 + k_avg - 2*q);

    % define mesh of m and rho values
    [m, rho] = meshgrid(-1:0.05:1, 0.001:0.05:0.501);

    % evaluate dm/dt and drho/dt vectors over mesh for current p and q
    dm = dm_dt(p, q, rho, 0.5.*(1+m), 0.5.*(1-m));
    drho = drho_dt(p, q, rho, 0.5.*(1+m), 0.5.*(1-m), k_avg);

    % calculate magnitudes of vectors
    magnitude = sqrt(dm.^2 + drho.^2);

    % normalize to unit vectors
    dm_unit = dm ./ magnitude; % unit vector in the m-direction
    drho_unit = drho ./ magnitude; % unit vector in the rho-direction

    % set desired length for all arrows
    arrow_length = 4; % can adjust
    dm_plot = dm_unit * arrow_length;
    drho_plot = drho_unit * arrow_length;

    % evaluate rho* for current p and q
    rho_star = ( (1 - p)*(k_avg - 2*q) - p ) / ( 2*(1-p)*(k_avg - q) );

    % update counter
    counter = counter + 1;

    % plot the arrows
    figure(2);
    set(gca,'fontsize',14)
    subplot(1,3,counter)
    hold on
    quiver(m, rho, dm_plot, drho_plot, 'LineWidth', 1.4, 'Color','#AF5F5A');
    xlabel('m');
    ylabel('\rho');

    xline(0, 'k', 'LineWidth', 2)
    plot(0,rho_star, 'k.', 'MarkerSize', 50)
    axis([-1 1 0 0.5])
    title(sprintf('p = pc, q = %.2f', q))

end

set(gca,'fontsize',14)
