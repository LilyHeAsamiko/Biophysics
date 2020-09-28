function Z = leftmost(U,Lambda)
    H = 0.5*log(2)+2.5*log(pi)-2.5*log(Lambda)-repmat(pi^2,6,1)./(2*Lambda)+0.5*Lambda;
    lU = log(U);
    X = exp(-pi^2/(2*Lambda));
    Z = X;
    K = Lambda/pi^2;
    % squeezing
    for t = 1:length(X)-1
        Z(t) = Z(t) - K(t)^(t^2-1);
        t = t+1;
        Z(t) = Z(t) + K(t)^(t^2-1);        
    end
    Z((reshape(H,6,1)+reshape(log(Z),6,1))<reshape(lU,6,1))= 0;
end