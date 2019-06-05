function E = solveE(X, L, T, beta, mu)

% to compute the noise tensor E

    mid = X - L + T/mu;
    Norder = ndims(mid);
    tempF = tenmat(mid,Norder);
    F = tempF.data;
    [row, col] = size(F);
    E = zeros(row, col);
    threshold = beta / mu;
    
    for i = 1:col
        t = norm(F(:,i));
        if t > threshold;
            E(:,i) = (t - threshold)/t * F(:,i);
        end
    end
    
    E = tenmat(E, tempF.rdims, tempF.cdims, tempF.tsize);
    E = tensor(E);
    E = E.data;

end