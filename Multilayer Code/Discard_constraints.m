function d=Discard_constraints(X,synf,syng,t,L_1,L_2)
% f constraints
d=0;
n=t*L_2;
for i=1:L_1
%    if sum(H(i,:))>0
        if VT(X((i-1)*n+1:i*n))~=synf(i)
            d=1;
            break;
        end
%    end
end
% g constraints
if d==0
    for j=1:L_2
%        if sum(H(:,j))>0
            if Check_Chunk(0,X,zeros(1,L_1),j,t,L_2)~=syng(j)
                d=1;
                break;
            end
%        end
    end
end