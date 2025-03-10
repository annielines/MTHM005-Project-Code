function[adj_matrix] = G_fixed_degree(n, d)

% function generating a random graph of n nodes each of fixed degree d:
% creates list of all nodes, each node appearing d times
% these nodes are then randomly paired (ensuring no self-loops or multi-edges)

% inputs:
% n = number of nodes
% d = degree of nodes (same for each node): must be even

% outputs:
% adj_matrix = adjacency matrix of generated network

% check inputted d is even and display error message if not
if mod(d,2) ~= 0

    disp('Change d so it is even.')

end

% initialise adjacency matrix with (n x n) zero matrix
adj_matrix = zeros(n);

% create list of nodes where each node appears d times
nodes = repelem(1:n, d);

% initialise counter
counter = 0;

% keep going until every node has been paired d times
while ~isempty(nodes)

    % initially set both nodes in pair as being the same (to enter while loop)
    node_1 = 1;
    node_2 = 1;

    % ensure pairing creates a valid edge (no self-loops or multi-edges)
    while (node_1 == node_2) || adj_matrix(node_1, node_2) == 1

        nodes = nodes(randperm(length(nodes))); % shuffle order of nodes

        % ensure last two nodes are not the same
        while nodes(end) == nodes(end-1)

            nodes = nodes(randperm(length(nodes))); % shuffle order of nodes

        end

        node_1 = nodes(1);
        node_2 = nodes(2);

        % add to counter counting number of pairings tried
        counter = counter + 1;

        if counter > n

            % return if remaining nodes are either all the same or already paired
            return

        end

    end

    % if valid, create the edge
    adj_matrix(node_1, node_2) = 1;
    adj_matrix(node_2, node_1) = 1;

    % remove the used nodes from the list
    nodes(1:2) = [];

    % reset counter
    counter = 0;
end

end