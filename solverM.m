
function M = solverM(M, U, Y, mu, gamma, p, R)
    % Update M

    numIter = 100;
    st = 0.005*max(abs(diag(M'*M)));  % sometimes a larger st will achieve better result
    D = ((M'*M + st*eye(R))^(1-p/2)); % advoid sigular value issue
    for i = 1:numIter
        M = (mu * U - Y) / (gamma * D + mu* eye(size(D))) ;
        D = ((M'*M + st*eye(R))^(1-p/2));
    end
end