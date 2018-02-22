%This function will fix all candidates using constraint f the primary 
%fixing, which needs received codeword, candidate for block boundaries
%and syndroms of VT f-checks, and also size of block and size of N
function [Y,blocks]=Prim_fix(X,can,syn,n)

[len,b]=size(can);
blocks=can;
N=b*n;
Y=zeros(len,N);

for i=1:len
    tt=length(X);
    Y(i,1:tt)=X;
    for j=1:b
        if can(i,j)==1
            blocks(i,j)=0;
            pointer=(j-1)*n+1-sum(blocks(i,1:j-1));
            pointer1=(j-1)*n+1-sum(can(i,1:j-1));

            S=VT_decoder(X(pointer1:pointer1+n-2),syn(j));

            Y(i,1:tt+1)=[Y(i,1:pointer-1),S,X(pointer1+n-1:end)];

            tt=tt+1;
        end
    end
end