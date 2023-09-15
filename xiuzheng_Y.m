%%%%�������Ӿ�����(��������������д���)
function Y = xiuzheng_Y(Y0,Su,Sv)  %arΪ˥��ϵ��

Yh0 = tenmat(Y0,1);  Yh = double(Yh0);
Yv0 = tenmat(Y0,2);  Yv = double(Yv0);
K = 50;
ar = 1; 

flagh = sum(Yh,2);    ind_h0 = find(flagh==0);    ind_h1 = find(flagh>0);   
flagv = sum(Yv,2);   ind_v0 = find(flagv==0);    ind_v1 = find(flagv>0);     



%%
Nh = size(Yh,1);
Su(1:(Nh+1):end) = 0;      %�Խ��߱�Ϊ0

Nv = size(Yv,1);
Sv(1:(Nv+1):end) = 0; 

Su(ind_h0,ind_h0) = 0;       %���������ʹ�л�����ƶ�Ϊ0
Sv(ind_v0,ind_v0) = 0;

[~,indh] = sort(Su,2,'descend');

ar = ar.^(0:K-1);
for i=1:length(ind_h0)    
    near_inh = indh(ind_h0(i),1:K);
    if sum((ar.*Su(ind_h0(i),near_inh)))>0
        Yh(ind_h0(i),:) = (ar.*Su(ind_h0(i),near_inh))*Yh(near_inh,:)/sum((ar.*Su(ind_h0(i),near_inh)));
    end
end

Yh =  tenmat(Yh,Yh0.rdims,Yh0.cdims,Yh0.tsize);
Yh = tensor(Yh);

[~,indv] = sort(Sv,2,'descend');
ar = ar.^(0:K-1);
for j=1:length(ind_v0)   
    near_inv = indv(ind_v0(j),1:K);
    if sum((ar.*Sv(ind_v0(j),near_inv)))>0
        Yv(ind_v0(j),:) = (ar.*Sv(ind_v0(j),near_inv))*Yv(near_inv,:)/sum((ar.*Sv(ind_v0(j),near_inv)));
    end
end 
Yv =  tenmat(Yv,Yv0.rdims,Yv0.cdims,Yv0.tsize);
Yv = tensor(Yv);

Y = double((Yh+Yv)/2);

end