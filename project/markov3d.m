function [h, E_steps, E] = markov3d(P,steps)
if steps < 25 
    b =2*steps+1;
else 
    b = 50;
end
a=zeros(b*b*b,3); 
for i = 1:b
    a((b*b*(i-1)+1):(b*b*(i-1)+b*b),1)=i;
    a((b*(i-1)+1):(b*(i-1)+b),2)=i;
end
a(:,2)=repmat(a(1:b*b,2),b,1);
z = 1:b;
a(:,3)=repmat(z',b*b,1);
s0 = b*b*steps+b*steps+steps; 
id = find(a(:,1)~=1);
id2 = find(a(id,1)~=b);
id2 = id(id2);
id3 = find(a(id2,2)~=1);
id3 = id2(id3);
id4 = find(a(id3,2)~=b);
id4 = id3(id4);
id5 = find(a(id4,3)~=1);
id5 = id4(id5);
id6 = find(a(id5,3)~=b);
id6 = id5(id6);
id7 = [];
for i = 1: size(P,1)
    id7 = find(a(id6,1)==P(i,1));
%    Id7 = [Id7; ID7];
    id7 = id6(id7);
    id8 = find(a(id7,2)==P(i,2));
    id8 = id7(id8);
    id9(i,:) = {find(a(id8,3)==P(i,3))};
    id10(i,:) = {id8(id9{i})};
    id7 = [];
    id8 = [];
end
id11 = [];
for i = 1:length(id10)
    id11 = [id11;id10{i}];
end
if id11
    id12 = [];
    id11 = unique(id11);
    for i = 1:length(id11)
        id12 = [id12;find(id6 ~= id11(i))];
    end
    id12 = id6(id12);
    si=find(id12==s0);
    s = 1:(length(id12)+1);
    ss=[id12(si:length(id12));id12(1:si-1)];
else
    si=find(id6==s0);
    s = 1:(length(id6)+1);
    ss=[id6(si:length(id6));id6(1:si-1)];
end

ind = sparse(length(s),length(s));
ind(length(s),length(s))=1;

for  i=1:(length(ss)-1)
    ii = find(a(ss,1)==(a(ss(i),1)-1));
    ii =ss(ii);
    ii2 = find(a(ii,2)==a(ss(i),2));
    ii2 = ii(ii2);
    ii3 = find(a(ii2,3)==a(ss(i),3));
    ii3 = ii2(ii3);
    iii = find(a(ss,1)==(a(ss(i),1)+1));
    iii =ss(iii);
    iii2 = find(a(iii,2)==a(ss(i),2));
    iii2 = iii(iii2);
    iii3 = find(a(iii2,3)==a(ss(i),3));
    iii3 = iii2(iii3);
    
    ii4 = find(a(ss,2)==(a(ss(i),2)-1));
    ii4 =ss(ii4);
    ii5 = find(a(ii4,1)==a(ss(i),1));
    ii5 = ii4(ii5);
    ii6 = find(a(ii5,3)==a(ss(i),3));
    ii6 = ii5(ii6);
    iii4 = find(a(ss,2)==(a(ss(i),2)+1));
    iii4 =ss(iii4);
    iii5 = find(a(iii4,1)==a(ss(i),1));
    iii5 = iii4(iii5);
    iii6 = find(a(iii5,3)==a(ss(i),3));
    iii6 = iii5(iii6);
    
    ii7 = find(a(ss,3)==(a(ss(i),3)-1));
    ii7 =ss(ii7);
    ii8 = a(ii7,1)==a(ss(i),1);
    ii8 = ii7(ii8);
    ii9 = find(a(ii8,2)==a(ss(i),2));
    ii9 = ii8(ii9);
    iii7 = find(a(ss,3)==(a(ss(i),3)+1));
    iii7 =ss(iii7);
    iii8 = find(a(iii7,1)==a(ss(i),1));
    iii8 = iii7(iii8);
    iii9 = find(a(iii8,2)==a(ss(i),2));
    iii9 = iii8(iii9);
end

idcord = [ii3(1) ii6(1) ii9(1) iii3(1) iii6(1) iii9(1)];
for i=1:length(idcord)/3
    ind(i,(idcord((i-1)*3+1):(idcord((i-1)*3+3))))=1/length(idcord)/3;    
end
%%

%for i = 11:16
i =15;
mc() = dtmc(ind((i-1)*1000+1:(i-1)*1000+1000,(i-1)*1000+1:(i-1)*1000+1000));
figure;
%subplot(2,3,i-10)
graphplot(mc(i-10));
%end

T = ind(1:length(ss),1:length(ss));
E =(speye(length(ss))-T)\speye(length(ss))*ones(length(ss),1)
%[AA, count] = blockmaker(T,e);
%% divide into small blocks and inverse
% d = length(T)/2;
% count = 1;
% while d>e
%     d = round(d/2);
%     count = count+1;
% end
% m = round(length(ss)/d);
% nn = length(ss)-d*m;
% if nn<0
%     for i =1:m-1
%          for j = 1:m-1
%              if i == j
%                  AA((i-1)*m+j,:)={1-T((i-1)*d+1:(i-1)*d+d,(j-1)*d+1:j*d)};
%                  AA((i-1)*m+m,:)={1-T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d+1:(j+1)*d+nn)};
%                  AA(m*(m-1)+j,:)={1-T((i-1)*d+d+1:(i+1)*d+nn,(j-1)*d+1:(j-1)*d+d)};
%              else
%                  AA((i-1)*m+j,:)={T((i-1)*d+1:(i-1)*d+d,(j-1)*d+1:j*d)};
%                  AA((i-1)*m+m,:)={T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d+1:(j+1)*d+nn)};
%                  AA(m*(m-1)+j,:)={T((i-1)*d+d+1:(i+1)*d+nn,(j-1)*d+1:(j-1)*d+d)};
%              end
%         end
%         if i == j    
%             AA(m*m,:)={1-T((m-2)*d+d+1:m*d+nn,(m-2)*d+d+1:m*d+nn)};
%         else
%             AA(m*m,:)={T((m-2)*d+d+1:m*d+nn,(m-2)*d+d+1:m*d+nn)};
%         end
%     end
% else
%     for i =1:m
%         for j = 1:m
%             if i == j
%                 AA((i-1)*m+j,:)={1-T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d:j*d)};
%                 AA((i-1)*m+m,:)={1-T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d+1:(j+1)*d+nn)};
%                 AA(m*(m-1)+j,:)={1-T((i-1)*d+d+1:(i+1)*d+nn,(j-1)*d+1:(j-1)*d+d)};
%             else
%                 AA((i-1)*m+j,:)={T((i-1)*d+1:(i-1)*d+d,(j-1)*d+1:j*d)};
%                 AA((i-1)*m+m,:)={T((i-1)*d+1:(i-1)*d+d,(j-1)*d+d+1:(j+1)*d+nn)};
%                 AA(m*(m-1)+j,:)={T((i-1)*d+d+1:(i+1)*d+nn,(j-1)*d+1:(j-1)*d+d)};
%             end
%         end
%         if i == j
%             AA(m*m,:)={1-T((m-2)*d+d+1:m*d+nn,(m-2)*d+d+1:m*d+nn)};
%         else
%             AA(m*m,:)={T((m-2)*d+d+1:m*d+nn,(m-2)*d+d+1:m*d+nn)};
%         end
%     end
% end

%TT = blockinverse(AA,m/2,m,count,1,e);         
%% times column vector ones
% for i = 1:length(ss)
%     E(i) = TT(i,:)*ones(length(ss),1);
% end

%%
%E = ones(length(ss),1)\(eye(length(ss))-T);
E_steps = max(E)
if E_steps > steps
    h = 0;
else
    
    h = 1;
end
end