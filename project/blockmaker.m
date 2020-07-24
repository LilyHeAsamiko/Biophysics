function [AA, count] = blockmaker(T,e)
%% divide into small blocks and inverse
d = length(T)/2;
count = 1;
while d>e
    d = round(d/2);
    count = count+1;
end
m = round(length(ss)/d);
nn = length(ss)-d*m;
if nn<0
    for i =1:m-1
         for j = 1:m-1
             if i == j
                 AA((i-1)*m+j,:)={1-T((i-1)*d+1:(i-1)*d+d,(j-1)*d+1:j*d)};
                 AA((i-1)*m+m,:)={1-T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d+1:(j+1)*d+nn)};
                 AA(m*(m-1)+j,:)={1-T((i-1)*d+d+1:(i+1)*d+nn,(j-1)*d+1:(j-1)*d+d)};
             else
                 AA((i-1)*m+j,:)={T((i-1)*d+1:(i-1)*d+d,(j-1)*d+1:j*d)};
                 AA((i-1)*m+m,:)={T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d+1:(j+1)*d+nn)};
                 AA(m*(m-1)+j,:)={T((i-1)*d+d+1:(i+1)*d+nn,(j-1)*d+1:(j-1)*d+d)};
             end
        end
        if i == j    
            AA(m*m,:)={1-T((m-2)*d+d+1:m*d+nn,(m-2)*d+d+1:m*d+nn)};
        else
            AA(m*m,:)={T((m-2)*d+d+1:m*d+nn,(m-2)*d+d+1:m*d+nn)};
        end
    end
else
    for i =1:m
        for j = 1:m
            if i == j
                AA((i-1)*m+j,:)={1-T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d:j*d)};
                AA((i-1)*m+m,:)={1-T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d+1:(j+1)*d+nn)};
                AA(m*(m-1)+j,:)={1-T((i-1)*d+d+1:(i+1)*d+nn,(j-1)*d+1:(j-1)*d+d)};
            else
                AA((i-1)*m+j,:)={T((i-1)*d+1:(i-1)*d+d,(j-1)*d+1:j*d)};
                AA((i-1)*m+m,:)={T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d+1:(j+1)*d+nn)};
                AA(m*(m-1)+j,:)={T((i-1)*d+d+1:(i+1)*d+nn,(j-1)*d+1:(j-1)*d+d)};
            end
        end
        if i == j
            AA(m*m,:)={1-T((m-2)*d+d+1:m*d+nn,(m-2)*d+d+1:m*d+nn)};
        else
            AA(m*m,:)={T((m-2)*d+d+1:m*d+nn,(m-2)*d+d+1:m*d+nn)};
        end
    end
end

end