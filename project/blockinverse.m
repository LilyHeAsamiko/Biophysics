function BB = blockinverse(AA,D,L,count,k,e)
%m = 2^count;
% AAA_=cell(4);
% AAA_(1,:)={AA{k}};% A_(k,:)
% AA_(2,:)={AA{k+D}};% A_(k+D,:)
% AA_(3,:)={AA{k+D*L}};% A_(k+D*L,:)
% AA_(4,:)={AA{k+D+D*L}};% A_(k+D+D*L,:)
% for i = 1:D-1
%     for j = 1:D-1
%         A_(1,:)={[A_{1},AA{k+(i-1)*L+j}]};
%         A_(2,:)={[A_{2},AA{k+(i-1)*L+j+D}]};
%         A_(3,:)={[A_{3},AA{k+D*L+(i-1)*L+j}]};
%         A_(4,:)={[A_{4},AA{k+D*L+(i-1)*L+j+D}]};
%     end
% end
A_=cell(L);
for i = 1:L
    A_(i,:)={AA{(i-1)*L+1}};
end
for i = 1:2*L
AA_(i,:)={A_{(i-1)*D+1}};% A_(k,:)
end
for i = 1:2+L
    for j = 1:D-1
        AA_(i,:)={[AA_{i},AA_{i+j}]};
    end
end

for i = 1:D-1
    AAA(1,:)={[AA_{1};AA_{1+2*i}]};
    AAA(2,:)={[AA_{2};AA_{2+2*i}]};
    AAA{3,:}={[AA_{2*D};AA_{2*D+2*i}]};
    AAA{4,:}={[AA_{2*D+1};AA_{2*D+1+2*i}]};
end

count = count-1
D = round(D/2)
if count>0
    size(AAA_{1})
    size(AAA_{2})
    size(AAA_{3})
    size(AAA_{4})
%    v = A_{4}-A_{3}*transpose(blockinverse(AA,D,L,count,k,e))*A_{2};
    v = AAA_{4}-(AAA_{3}.*blockinverse(AA,D,L,count,k,e)).*AAA_{2};
%    [vv, count] = blockmaker(v,e);
%    inv_v=blockinverse(v,D,L,count,k);
    
    B1 = blockinverse(AAA,D,L,count,k,e)+blockinverse(AAA,D,L,count,k,e)*AAA_{2}*blockinverse(v,D,L,count,k,e)*AAA_{3}*blockinverse(AA,D,L,count,k,e);
    B2 = -blockinverse(AA,D,L,count,k,e)*AAA_{2}*blockinverse(v,D,L,count,k,e);
    B3 = -blockinverse(v,D,L,count,k,e)*AAA_{3}*blockinverse(AA,D,L,count,k,e);
    B4 = blockinverse(v,D,L,count,k,e);
    diplay("B:")
    size(B1)
    size(B2)
    size(B3)
    size(B4)
else
    display('last:')
    size(AAA_{1})
    size(AAA_{2})
    size(AAA_{3})
    size(AAA_{4})
    AAA=[AAA_{1},AAA_{2};AAA_{3},AAA_{4}];
    AAA(1,:) = {AAA(1:size(AAA_{1},1),1:size(AAA_{1},1))};
    AAA(2,:) = {AAA(1:size(AAA_{1},1),(size(AAA_{1},1)+1):size(AAA,2))};
    AAA(3,:) = {AAA((size(AAA_{1},1)+1):size(AAA,1),1:(size(AAA,2)-size(AAA_{1},1)))};
    AAA(4,:) = {AAA((size(AAA_{1},1)+1):size(AAA,1),(size(AAA,2)-size(AAA_{1},1)+1):size(AAA_{1},2))};
    v1 = AAA_{1}-AAA_{2}*(AAA_{3})'*inv(AAA_{4});
    inv1 =inv(v1);
    B1 = inv1;
    B2 = -inv1*inv(AAA_{4})*AAA_{2};
    B3 = -inv(AAA_{4})*inv1*AAA_{3};
    B4 = inv(AAA_{4})+inv(AAA_{4})*inv(AAA_{4})*inv1*(AAA_{2})*AAA_{3}';
end
BB = [B1,B2;B3,B4];
end
