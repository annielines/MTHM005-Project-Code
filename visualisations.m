%% typical network configurations at each stable steady state
clear

N = 100; % set number of nodes
K = 4; % set degree

% set initial condition m = 0
initial_1 = 0.5;

% set maximum number of steps of model until it terminates
max_real = 20000;

% generate a graph of fixed degree K with N nodes
g = G_fixed_degree(N, K);

%% coexistence: take (p,q) = (0.2, 0.5)

% set parameter values
p = 0.2;
q = 0.5;

% implement coevolutionary nonlinear voter model
[coex_adj, ~, final_0, final_1] = coev_nonlinear_voter_model(g, p, q, initial_1, max_real);

% convert output into a graph
coex_g = graph(coex_adj);

% set blue nodes for opinion state 0, red nodes for opinion state 1
node_colors(final_0, :) = repmat([0, 0, 1], length(final_0), 1);
node_colors(final_1, :) = repmat([1, 0, 0], length(final_1), 1);

% plot the coexistence graph with coloured nodes
figure(1)
coex_plot = plot(coex_g, 'Layout', 'force');
coex_plot.NodeColor = node_colors;

%% single component global consensus: take (p,q) = (0.2, 2)

% set parameter values
p = 0.2;
q = 2;

% implement coevolutionary nonlinear voter model
[con_adj, ~, final_0, final_1] = coev_nonlinear_voter_model(g, p, q, initial_1, max_real);

% convert output into a graph
con_g = graph(con_adj);

% set blue nodes for opinion state 0, red nodes for opinion state 1
node_colors(final_0, :) = repmat([0, 0, 1], length(final_0), 1); % blue for 0
node_colors(final_1, :) = repmat([1, 0, 0], length(final_1), 1); % red for 1

% plot the coexistence graph with coloured nodes
figure;
con_plot = plot(con_g, 'Layout', 'force');
con_plot.NodeColor = node_colors;

%% fragmentated consensus,: take (p,q) = (0.8, 0.5)

% set parameter values
p = 0.8;
q = 0.5;

% implement coevolutionary nonlinear voter model
[frag_adj, ~, final_0, final_1] = coev_nonlinear_voter_model(g, p, q, initial_1, max_real);

% convert output into a graph
frag_g = graph(frag_adj);

% set blue nodes for opinion state 0, red nodes for opinion state 1
node_colors(final_0, :) = repmat([0, 0, 1], length(final_0), 1); % blue for 0
node_colors(final_1, :) = repmat([1, 0, 0], length(final_1), 1); % red for 1

% plot the fragmented graph with coloured nodes
figure;
frag_plot = plot(frag_g, 'Layout', 'force');
frag_plot.NodeColor = node_colors;
