function[adj_matrix, consensus_time, opinions_0, opinions_1, conflicting_edges, rho_time_series,  opinion_1_time_series] = coev_nonlinear_voter_model(adj_matrix, p, q, initial_1, max_steps)

% function for coevolving nonlinear voter model (binary discrete opinions)

% here, indicate different opinion states with either 1 or 0

% inputs:
% adj_initial = adjacency matrix for initial social network
% 1_initial = initial proportion of nodes with opinion 1
% p = probability of re-wiring
% q = degree of nonlinearity
% max_steps = maximum number of update steps before terminating model

% outputs:
% adj_matrix = adjacency matrix for final state of social network
% consensus_time = number of update steps taken to reach consensus
% opinions_0 = vector of indexes of all nodes with final opinion 0
% opinions_1 = vector of indexes of all nodes with final opinion 1
% conflicting_edges = matrix indicating nodes that are connected by an active/conflicting edge
% rho_time_series = vector of rho at each update step
% opinion_1_time_series = vector of proportion of nodes with opinion 1 at each update step

% find number of nodes in social network
N = length(adj_matrix);

% find number of edges in social network
edges = 0.5*(sum(sum(adj_matrix)));

% calculate the desired number of nodes to be assigned opinion 1
num_ones = round(initial_1 * N);

% the rest of the nodes will have opinion 0
num_zeros = N - num_ones;

% create a vector containing 'num_ones' 1s and 'num_zeros' 0s
opinions = [ones(1, num_ones), zeros(1, num_zeros)];

% shuffle the vector to randomise opinion assignement
opinions = opinions(randperm(N));

% find the nodes that are initially assigned each opinion
opinions_0 = find(opinions == 0); % opinion 0
opinions_1 = find(opinions == 1); % opinion 1

% initialise counter tracking time to reach consensus
consensus_time = 0;

% initialise opinion difference matrix (*)
D = abs(opinions - transpose(opinions));

% multiply the adjacency matrix with D, element-wise (**)
conflicting_edges = D.*adj_matrix;

% track proportion of each vote over time
opinion_0_time_series(1) = num_zeros / N; % initialise for 0s over time
opinion_1_time_series(1) = num_ones / N; % initialise for 1s over time

% track rho over time
rho_time_series(1) = (sum(sum(conflicting_edges)) / 2) / edges;

% continue going until a consensus has been reached
while sum(sum(conflicting_edges)) ~= 0

    % select a node at random
    node_i = randi(N);

    % ensure this node is connected to at least one other node
    while sum(adj_matrix(node_i, :)) == 0
        node_i = randi(N);
    end

    % count number of active edges adjacent to node i
    a_i = sum(conflicting_edges(:, node_i));

    % count total degree of node i
    k_i = sum(adj_matrix(:, node_i));

    % calculate pi for current node
    pi_i = a_i / k_i;

    % with probability pi_i...
    if rand < pi_i^q

        % find all active edges
        active_edges = find(conflicting_edges(node_i, :) == 1);

        % select a random neighbouring node connected by an active edge
        j = randi([1 length(active_edges)]);
        node_j = active_edges(j);

        % rewire with probability p
        if rand < p

            % i.e. form new connection with node k with same opinion
            % (if node i has a unique opinion, or is already connected to
            % every other node with same opinion, skip rewiring step)

            opinions_0_temp = find(opinions == 0); % find all nodes with opinion 0
            opinions_1_temp = find(opinions == 1); % find all nodes with opinion 1

            % remove node i from these lists
            opinions_0_temp(opinions_0_temp == node_i) = [];
            opinions_1_temp(opinions_1_temp == node_i) = [];

            % if node i has opinion 0 and is not the last opinion 0
            if (opinions(node_i) == 0) && (length(opinions_0_temp) ~= 0)

                % select random node k with opinion 0
                node_k_index = randi([1 length(opinions_0_temp)]);
                node_k = opinions_0_temp(node_k_index);

                % ensure node i and node k not already connected
                while (adj_matrix(node_i, node_k) == 1) && (length(opinions_0_temp) ~= 0)

                % remove node k from potential new connections
                opinions_0_temp(opinions_0_temp == node_k) = [];

                if (length(opinions_0_temp) ~= 0)

                    % select new random node k with opinion 0
                    node_k_index = randi([1 length(opinions_0_temp)]);
                    node_k = opinions_0_temp(node_k_index);

                end

                end

                if(length(opinions_0_temp) ~= 0)

                    % delete connection between node i and node j
                    adj_matrix(node_i, node_j) = 0;
                    adj_matrix(node_j, node_i) = 0;

                    % form new connection between node i and node k
                    adj_matrix(node_i, node_k) = 1;
                    adj_matrix(node_k, node_i) = 1;

                end

            % if node i has opinion 1 and is not the last opinion 1
            elseif (opinions(node_i) == 1) && (length(opinions_1_temp) ~= 0)

                % select random node k with opinion 1
                node_k_index = randi([1 length(opinions_1_temp)]);
                node_k = opinions_1_temp(node_k_index);

                % ensure node i and node k not already connected
                while (adj_matrix(node_i, node_k) == 1) && (length(opinions_1_temp) ~= 0)

                % remove node k from potential new connections
                opinions_1_temp(opinions_1_temp == node_k) = [];

                if (length(opinions_1_temp) ~= 0)

                    % select new random node k with opinion 1
                    node_k_index = randi([1 length(opinions_1_temp)]);
                    node_k = opinions_1_temp(node_k_index);

                end

                end

                if(length(opinions_1_temp) ~= 0)

                    % delete connection between node i and node j
                    adj_matrix(node_i, node_j) = 0;
                    adj_matrix(node_j, node_i) = 0;

                    % form new connection between node i and node k
                    adj_matrix(node_i, node_k) = 1;
                    adj_matrix(node_k, node_i) = 1;

                end

            end


        % copy opinion of node j with probability (1-p)
        else

            opinions(node_i) = opinions(node_j);

        end

    end

    % update counter
    consensus_time = consensus_time + 1;

    % update opinion difference matrix (*)
    D = abs(opinions - transpose(opinions));

    % multiply the adjacency matrix with D, element-wise (**)
    conflicting_edges = D.*adj_matrix;

    % find the nodes that now have each opinion
    opinions_0 = find(opinions == 0); % opinion 0
    opinions_1 = find(opinions == 1); % opinion 1

    % update tracker for proportion of each opinion over time
    opinion_0_time_series(consensus_time) = length(opinions_0) / N;
    opinion_1_time_series(consensus_time) = length(opinions_1) / N;

    % update tracker for rho over time
    rho_time_series(consensus_time) = (sum(sum(conflicting_edges)) / 2) / edges;

    if consensus_time > max_steps

        break

    end

end

% notes:

% (*) (i,j)th element is a 1 if node i and node j have different opinion, else 0

% (**) only elements that are 1s in both adjacency matrix (i.e. connected) and
% in D (i.e. different opinions) are 1s in this product

end