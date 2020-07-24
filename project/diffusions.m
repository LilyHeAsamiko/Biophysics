function [ss,h] = diffusions(N,Steps,dx)
    if Steps < 25 
        b =2*Steps+1;
    else 
        b = 50;
    end
    P = randi(b,N,3);
    B = [P-[1,0,0];P+[1,0,0];P-[0,1,0];P+[0,1,0];P-[0,0,1];P+[0,0,1]];
    [ss,h] = markov3d2(B, Steps, dx)
end

    