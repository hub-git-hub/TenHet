function [L, errL] = TenHet(X, Omega, A, beta, gamma, lambda, p, R, maxIter, tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This code solves the tensor completion with mutiple heterogeneous 
% data source via Alternation Direction Method of Multipliers (ADMM)
% 
% INPUT
% X: observed 3-way tensor
% Omega: observed indexes of tensor
% A: mutiple side information for mode of tensor
%  A{i}{j} denotes the j-th view information for mode i

% beta: regularized parameter for noise tensor
% gamma: regularized parameter for Schatten p-norm of latent factors
% lambda: regularized parameter for coupled tensor-matrix decomposition
% R: the rank of tensor
% maxIter, tol: user-defined max number of iterations and tolerance

%
% OUTPUT
% L:  full recovery tensor of X
% errL: a list of error reflect the convergence of L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default parameters
if nargin < 10
    tol = 1e-5;  
end

if nargin < 9
    maxIter = 50; 
end


if nargin < 8
    R = 20;
end

if nargin < 7
    p = 0.5;
end


mu = 1e-3;

% initialize the complete tensor
L = X;

% initialize the error tensor E
E = 0*rand(size(X));

% initialize error list to check convergence of L
errL = zeros(maxIter, 1);

% number of dimensions of the observed tensor
ndims_X = ndims(X);

% norm of the observed tensor
normX   = norm(X(:));

% size of the observed tensor
X_size  = size(X);

% number of view;
ni = length(A{1});


% initialize latent matrices {U1, U2, ..., Un} with random number
% initialize auxiliary matrices {M1, M2, ..., Mn} as {U1, U2, ..., Un}
% initialize Lagrange multipler {Y1, Y2, ..., Yn} with zeros
size_II = []; 
for i = 1:ndims_X       
    U{i} = rand(X_size(i),  R);
    M{i} = U{i};
    Y{i} = zeros(X_size(i), R);
    size_II = [size_II, R];
end

index_temp = ones(1, R);
II = tendiag(index_temp, size_II);


% initialize latent matrices G for side informaiton A,B,C
% initialize auxiliary variable Q = G; 
% initialize lagrangian multipler Z = 0  and scale matrix S;
for i = 1:ndims_X
    for j = 1:ni
        dimm = size(A{i}{1});
        G{i}{j} = rand(dimm(1),R);
        Q{i}{j} = G{i}{j};
        Z{i}{j} = zeros(dimm(1),R);
        S{i}{j} = diag(sum(G{i}{j}));
    end
end

% initialize lagrangian multipler T
 T = zeros(X_size);
    
% record the previous estimated results
L_pre = X;





% % --------------------------- Iteration Scheme ----------------------------
for iter = 1: maxIter   
    
    % update M{i}
    for i = 1:ndims_X  
        M{i} = solverM(M{i}, U{i}, Y{i}, mu, gamma, p, R);
    end

    
    % update U{i}
    for i = 1:ndims_X    
        midT = tensor(II);
        % calculate Khatric-Rao product of U(1), ..., U(i-1),U(i+1), ...,U(n)
        for m = 1:ndims_X
            if m == i
            continue;
            end
            midT = ttm(midT, U{m}, m);
        end
        unfoldD_temp = tenmat(midT, i);  
        temp_B = 2*unfoldD_temp.data*unfoldD_temp.data';
        temp_B = temp_B + ( 2 * lambda * ni + mu)*eye(R,R);
        temp_B = temp_B + 0.00001*eye(size(temp_B));
        temp_C = tenmat(L, i);
        temp_D = 2 * temp_C.data*unfoldD_temp.data';

        for k = 1:ni

            temp_D = temp_D + 2*lambda*G{i}{k}*S{i}{k};
        end
        temp_D = temp_D + mu * M{i} + Y{i};

        U{i}   = temp_D /temp_B;
    end
    clear unfoldD_temp temp_B temp_Z temp_C  ;


    % update G{i}{j},  Q{i}{j} and Z{i}{j}
    for i = 1:ndims_X
        for j = 1:ni
            temp_A = 2 * A{i}{j}' * Q{i}{j};
            temp_A = temp_A + 2*U{i}* S{i}{j}';
            temp_A = temp_A + mu *Q{i}{j} + Z{i}{j};

            temp_B = 2*Q{i}{j}'*Q{i}{j} +  2*S{i}{j}*S{i}{j}' + mu * eye(R,R);
            temp_B = temp_B + 0.00001*eye(size(temp_B));
            G{i}{j} = temp_A / temp_B;
            S{i}{j} = diag(sum(G{i}{j}));

            temp_A = 2 *A{i}{j}*G{i}{j} + mu * G{i}{j} - Z{i}{j};
            temp_B = 2 *G{i}{j}'*G{i}{j} + mu * eye(R,R);
            temp_B = temp_B + 0.0001*eye(size(temp_B));
            Q{i}{j} = temp_A / temp_B;

            Z{i}{j} = Z{i}{j} + mu*(Q{i}{j} -G{i}{j});
        end
    end
    clear temp_A  temp_B;


    % update  Y{i}
        for i = 1:ndims_X  
            Y{i} = Y{i} + mu * (M{i} - U{i});
        end


    % update L
    midT = tensor(II);
    midT = ttm(midT, U, [1:ndims_X]);  
    L = midT.data;
    L = (2*L + mu*X - mu*E + T)/(mu+2);
    
    
    
    % update the noise tensor E
    E = solveE(X, L, T, beta, mu);

    % whos
    
    % update T
    T = T + mu*(X - L - E);


    % checking the stop criteria
    stopC = norm(L_pre(:) - L(:))/norm(L_pre(:));
    L_pre = L; 
    errL(iter) = stopC;
    L(Omega) = X(Omega);

    if stopC < tol
       fprintf('the difference less the tolerance, and it stops\n');
       break;
    end  
end
end

