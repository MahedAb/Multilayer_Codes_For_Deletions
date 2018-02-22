function [X,S]=ready_Lin(X,H,t) 
% n is block length and t is chunk length
%deletions will replaced by -1 (a way to show erasures)
%S will save the string for that chunk to check it afterwards with erasures
%solutions, third element of each row shows length of original string,
%excess elements of S are zeros.

[L_1,L_2]=size(H);
n=t*L_2;
S=zeros(nnz(H),t+3);
%excess=sum(sum(H));
%X=X(1:end-excess);
l=0;
for i=1:L_1
    if sum(H(i,:))>0
        for j=1:L_2
            if H(i,j)>0
                l=l+1;
                ind=(i-1)*n+(j-1)*t;
                S(l,1:3)=[i,j,(t-H(i,j))];
                S(l,4:(t-H(i,j))+3)=X(ind+1:ind+(t-H(i,j)));
                X=[X(1:ind),-1*ones(1,t),X(ind+(t-H(i,j))+1:end)];
            end
        end
    end
end