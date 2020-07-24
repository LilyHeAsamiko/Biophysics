function [i, h] = markov3d2(B, Steps, dx)
x = [26, 26, 26];
direction = [1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1];
for i = 1:Steps
    x = x+dx*direction(randi(6,1,1),1);
    if ((x(1) ==1) || (x(1) == 51) || (x(2) ==1) || (x(2) ==51) || (x(3) == 1) || (x(3) == 51))
        h = 0;
        break;
    else
        h=2;
    end
    for j = 1:length(B)
        if x==B(j,:)
           h = 1;
           break;
        end
    end
    if h==1
        break;
    end
end
end