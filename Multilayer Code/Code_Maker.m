function [X,fsyn,gsyn,del_bit,Y,L,r,A]=Code_Maker(t,L_1,L_2,k,lin)
n=t*L_2;
fsyn=floor((n+1)*rand(1,L_1));
A=rand(L_1,n)>0.5;
N=n*L_1;

X=zeros(1,N);
for i=1:L_1
    A(i,:)=VT_decoder(A(i,1:n-1),fsyn(i));
    X((i-1)*n+1:i*n)=A(i,:);
end

s=zeros(t*L_1,1);
gsyn=zeros(L_2,1);
for i=1:L_2
    for j=1:L_1
        s((j-1)*t+1:j*t)=A(j,(i-1)*t+1:i*t);
    end
    gsyn(i)=VT(s);
end

Y=X;
R=randperm(n*L_1);
R=sort(R(1:k));
del_bit=zeros(1,k);
L=R;
for i=1:k
    del_bit(i)=ceil(L(i)/n);
    if R(i)<length(X)
        X=[X(1:R(i)-1),X(R(i)+1:end)];
    else
        X=X(1:R(i)-1);
    end
    if i<k
        R(i+1:end)=R(i+1:end)-1;
    end
end

x=zeros(1,L_1);
for i=1:length(del_bit)
    x(del_bit(i))=x(del_bit(i))+1;
end
del_bit=x;

A=rand(lin,N)>0.5;
r=mod(A*Y',2);
