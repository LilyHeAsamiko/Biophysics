function R = ac(A)
%             R = zeros(length(A));
%             M = A'*A;
%             m = mean(mean(M));
%             N = M - m;
%             M2 = M - repmat(mean(M,1),length(A),1);
%             for i = 1:length(A)
%                 for j = [i+1:length(A),1:i]
% %                for j = i:length(A)
%                     R(i,j)= sum(M2(:,i).*M2(:,j))/std(M(:,i))*std(M(:,j));
% %                    R(i,j)= mean(cov(M(:,i),M(:,j)))/std(M(:,i))*std(M(:,j));
%                 end
%             end
%             R = R/repmat(mean(R,1),length(A),1);
    avgA = mean(A);
    c0 = cov(A);
    R = zeros(length(A),1);
    for k = 0:length(A)-1
        ck = 0;
        for i = 1:(length(A)-k)
            ck =ck + (A(i)-avgA)*(A(i+k)-avgA);
        end
        R(k+1) = ck/(length(A)-1)/c0;    
    end
end
        