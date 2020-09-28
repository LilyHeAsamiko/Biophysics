function Z = rightmost(U,lambda)
    X = exp(-0.5*lambda);
    Z = X;
    % squeezing
    for t = 1:length(X)-1
            Z(t) = Z(t) - (t+1)^2*X(t)^(t+1)^2-1;
            t = t+1;
            Z(t) = Z(t) + (t+1)^2*X(t)^(t+1)^2-1;        
    end
    Z(Z<U)= 0;
end