a = [25:49, 24:-1:1];
ind = zeros(length(a));
ind(length(a),length(a))=1;
for  i=1:length(s)-1
    if (i == find(a==49)) 
        ind(i,length(a))=0.5;
        ind(i,i-1)=0.5;
    end
    ind(i,i+1)=0.5;
    ind(i,find(a==(a(i)-1)))=0.5;
end

m_ind = zeros(length(a));
m_ind(:,length(a))=1;

mu_ind = zeros(length(a));
p = zeros(1,length(a));
for i = length(a):-1:1
    s = 0;
    for j = length(a):-1:i
        s = s+ind(j,j)
    end
    p(50-i)=sum(sum(ind(:,length(a):-1:i),1))-s
end
for i = 1:length(a)
    mu_ind(:,i)=p(i)/sum(p(:));
end



