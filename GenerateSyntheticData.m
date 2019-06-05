function [X, A,index] = GenerateSyntheticData(tenSize, K, n_views)

% generate A 3-way block diagonal tensor X  and 
% multiple block diagonal matrices A for each mode of tensor

% Input:
%	tenSize: the size of tensor
%	K: the number of clusters locating on the digonal of tensor
%	n_views: number of views of side information for each mode of tensor

% Output:
%	X: tensor e.g., tenSize X tenSize X tenSize
%	A: the multi-view side information, e.g. A{i}{j} denotes the j-th
%      view side information for mode-i of tensor X
%   index: the index of samples from block tensor
%          it use for spliting training and test data later

tenOrder = 3;
gap = tenSize / K;
Xdim = tenSize * ones(1, tenOrder);

X = rand(Xdim);  % noise tensor
index = rand(Xdim); % index for sample traning and test data later

for i = 1:tenOrder
	for j = 1:n_views
		A{i}{j} = 0.2 * i * rand(tenSize ,tenSize); % different level noise matrix 
	end
end


% Generate the block dignoal tensor X and  block diagonal matrix A with noise
diff = [1:gap:tenSize];
for a = 1: K
	for i = diff(a): (diff(a) + gap - 1)
		for j = diff(a): (diff(a) + gap - 1)
			for k = diff(a): (diff(a) + gap - 1)
				X(i,j,k) = 1 + rand(1, 1);    % emphasize the diagonal
				index(i,j,k) = 1- 0.1 * rand(1,1); % emphasize the diagonal
			end

		    for mode = 1: tenOrder
			    for view =1: n_views
				    A{mode}{view}(i,j) = 1- 0.5 * view * rand(1,1); % emphasize the diagonal
			    end
		    end
		end
	end
end


% make A symmetric matrices

for mode = 1: tenOrder
	for view =1: n_views
		A{mode}{view} = (A{mode}{view}' + A{mode}{view} ) / 2;
	end
 end

end


