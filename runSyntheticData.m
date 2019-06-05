clc;clear all;

addpath('tensor_toolbox')

tenSize = 100;
numcluster = 5;
n_views = 3;

%%% generate a 100X100X100 tensor X with 5 cluster on the block diagonal
%%% generate three views similarity matrices A for each mode of tensor X

[data, A, index] = GenerateSyntheticData(tenSize, numcluster, n_views);


% one sample of side information A{1}{1}
aa = figure(1)
subplot(1,3,1)
imshow(A{1}{1},[]); colormap jet; colorbar;
title('Low noise','FontSize',16)
subplot(1,3,2)
imshow(A{1}{2},[]); colormap jet; colorbar;
title('Medium noise','FontSize',16)
subplot(1,3,3)
imshow(A{1}{3},[]); colormap jet; colorbar;
title('High noise','FontSize',16)



%%%%----------------Main--------------------------
%%% parameters setting up
fraction = 0.5;                % fraction of missing data
Omega = index > fraction;     % observed index for training set

X = data;
X(logical(1-double(Omega))) = 0;  % observed tensor


%%% compete the observed tensor X
tic
  [L, errX_hat] = TenHet(X, Omega, A, 0.001, 0.01, 0.01, 0.5);
time = toc;

%%% calculate relative error
relErr = norm(L(:) - data(:), 'fro')/ norm(data(:), 'fro');
fprintf('The relative error of tensor L is: %f.\n', relErr);


%%% Convergence check
bb = figure(2)
plot(errX_hat(:))
xlabel('Number of iterations')
ylabel('Relative error')
title('Convergence')




saveas(aa,'Noises.png')
saveas(bb,'Convergence.png')


