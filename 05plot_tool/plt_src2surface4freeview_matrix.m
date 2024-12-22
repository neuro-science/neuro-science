function M = plt_src2surface4freeview_matrix(anatomical, node_ids, method, depth, alpha)
% % % written 22/08/24 by wp: generate weigt matrix for freeview to smooth data

	if nargin < 2
        error('We need the faces(tris)/vertices and source ids anyway!');
    end

    if nargin < 3 || isempty(method)
        method = 'Barycentric';
    end

    if strcmpi(method(1), 'b')	%	Barycentric
			% % % In this case, anatomical is Nx3 faces
		 
        if nargin < 4 || isempty(depth)
            depth = 21;	%propagands to how many steps
        end

        if nargin < 5 || isempty(alpha)
            alpha = 1;	%decay alpha / (step+1)
        end

        num_vertices = max(anatomical(:));
        num_data_nodes = length(node_ids);

        if num_vertices < 0 || any(node_ids > num_vertices)
            error('Size mismatch: check triangulation and node_ids.');
        end

        M = zeros(num_vertices, num_data_nodes);

        avg_neighbors = 7; % Estimated average number of neighbors per vertex

        for i = 1:num_data_nodes
            current_node = node_ids(i);
            visited = false(num_vertices, 1);
            visited(current_node) = true;
            queue = [current_node];
            steps = 0;

            while steps < depth && ~isempty(queue)
                steps = steps + 1;

                % Preallocate next_queue based on estimated size
                next_queue = zeros(length(queue) * avg_neighbors, 1);
                next_queue_idx = 1;

                % Precompute connected faces to optimize loop performance
                connected_faces = find(any(ismember(anatomical, queue), 2)); 

                for j = 1:length(queue)
                    current_vertex = queue(j);

                    % Process only the faces connected to the current vertex
                    relevant_faces = connected_faces(any(anatomical(connected_faces, :) == current_vertex, 2));
                    
                    for face_idx = relevant_faces'
                        face_vertices = anatomical(face_idx, :);
                        neighbors = face_vertices(face_vertices ~= current_vertex);

                        % Loop through each neighbor individually
                        for k = 1:length(neighbors)
                            neighbor = neighbors(k); % Get the individual neighbor

                            if ~visited(neighbor)
                                visited(neighbor) = true;
                                next_queue(next_queue_idx) = neighbor;
                                next_queue_idx = next_queue_idx + 1;

                                % Calculate the weight based on the distance using the decay factor alpha
                                weight = 1 / (steps + 1)^alpha;

                                % Assign the weight to the matrix
                                M(neighbor, i) = M(neighbor, i) + weight;
                            end
                        end
                    end
                end

                % Trim the unused part of the preallocated next_queue
                next_queue = next_queue(1:next_queue_idx - 1);

                queue = next_queue;
            end

            % Ensure the data point itself has full weight
            M(current_node, i) = 1;
        end
	 elseif strcmpi(method(1), 'g')	%	Gaussian
		% % % In this case, anatomical is Nx3 vertices
		 
		if nargin < 4 || isempty(depth)
			depth = 1000;	% cut-off to zero at this distance
		end

		if nargin < 5 || isempty(alpha)
			alpha = 0.02;	% reduce to 50% at this distance
		end
		 
		% % % compute the transform matrix	
		alpha2 = (alpha / sqrt(log(2))) .^ 2;
		dist2 = sum(bsxfun(@minus, permute(anatomical, [1 3 2]), ...
			permute(anatomical(node_ids, :), [3 1 2])).^2, 3);
		M = exp(-dist2 / alpha2);
		M(dist2 > depth.^2) = 0; % Apply cutoff
		M(isnan(M)) = 0; % Remove nan when divided by zeros
	
    else
        error('Unsupported interpolation method.');
    end
end
