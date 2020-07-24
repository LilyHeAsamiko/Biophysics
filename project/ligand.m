N = 1;
for n = 1:N
    ss = zeros(1,N);
    h = zeros(1,N);
    steps = zeros(1,N);
    s1 = zeros(1,N);
    NN = 1000;
    for j = 1:NN
        R = 100*10^(-12);
        A = 2.414*10^(-5);
        B = 247.8;
        C = 140;
        T = 37+273.15;
        eta = A*10^(B/(T-C));
        t = 1;
        kb = 1.38*10^(-23);
        D = kb*T/(6*pi*eta*R);
        dx = 1*10^(-9);
        dt = dx^2/(4*D);
        Steps = round(t/dt)+1;
        [ss(j,n), h(j,n)] = diffusions(N,Steps,dx)
        if h(j,n) == 1
            steps(j,n) = ss(j,n);
        elseif h(j,n) == 2
            s1(j,n) = ss(j,n);
            end44444444444
    end
end
p = sum(h==1,1)/1 % the nth term means for n enzymes, the possibility for one to be found in total NN times  
N0 = sum(h==0,1) %the nth term means for n enzymes, the ligand hit the edge N0 times 
N1 = sum(h==1,1) % the nth term means for n enzymes, one is found N1 times 
steps = sum(ss,1)/N1;% the nth term means for n enzymes, average successful steps
s1 = sum(s1,1)/N0;% the nth term means for n enzymes, average edge hitting steps
[n_N1, max_N1] = max(N1);
[n_N0, min_N0] = min(N0);
[n_steps, min_steps] = min(steps);
[n_s1,min_s1] = min(s1);
figure,
subplot(1,3,1)
plot(1:NN, p);
title('probability to find enzymes');
subplot(1,3,2)
plot(1:NN, N0);
title('times hitting the edge ');
subplot(1,3,3)
plot(1:NN,steps);
title('average steps managing to find enzymes');
[X,Y]=meshgrid(steps,N1);
figure,
pcolor(X,Y,p);
title('average steps and its possiblity to find enzymes');
