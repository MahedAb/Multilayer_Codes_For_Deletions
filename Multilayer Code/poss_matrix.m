%This function creat all possiblities for deletions in chunks based on
%deletions that left
%ind show that whether the VT constraint is satisfied or not
function A=poss_matrix(V,ind,ch)
t=length(V);
%t here is number of blocks thathave deletion left (more than one)

V=min(ch,V);
A=(0:V(1))';
for i=2:t
    [u,~]=size(A);
    A=repmat(A,V(i)+1,1);
    for j=0:V(i)
        A(j*u+1:(j+1)*u,i)=j*ones(u,1);
    end
end
if ind==0
    A=A(2:end,:);
elseif ind==1
    if V(1)>0
    A=[A(1,:);A(3:end,:)];
    end
    s=V(1)+1;
    j=1;
    for i=2:t
        j=j+1-(V(i-1)==0);
        if V(i)>0
            A=[A(1:(s-(j-1)),:);A(s+2-(j-1):end,:)];
            s=s*(V(i)+1);
        end
    end
end