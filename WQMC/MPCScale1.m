function [W,Z,Ca,PhiA,GammaA] = MPCScale1(A,B,C,Np)
D = 0;
A_d = A;
B_d = B;
C_d = C;
u = size(B_d,2);
n =size(A_d,1);
p =size(C_d,1);
% build PhiA;
PhiA = [A_d sparse(n,p);
    C_d*A_d speye(p)];
% build GammaA;
GammaA = [B_d;C_d*B_d];
PhiA = sparse(PhiA);
GammaA = sparse(GammaA);
Ca = [sparse(p,n) speye(p)];
[mCa,nCa] = size(Ca);
%build W;
W_Transpose = sparse(nCa,mCa*(Np+1));
temp = Ca';
PhiA_T = PhiA';
for i = 1:Np+1
    W_Transpose(:,(mCa*(i-1)+1):mCa*i) =  temp;
    temp = PhiA_T*temp;
end

W = W_Transpose(:,(mCa+1):end)';

Z = sparse(mCa*Np,u*Np);

temp = W_Transpose(:,1:mCa*Np);
temp = temp' * GammaA;
 Z(1:mCa*Np,1:u) = temp;
for i = 2:Np
    Z((i-1)*mCa +1:mCa*Np,(i-1)*u + 1: i*u) = temp(1:(Np-i+1)*mCa,:);
end

end
