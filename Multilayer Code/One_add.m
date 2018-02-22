function Pos=One_add(X)
n=length(X);
Pos=zeros(n+2,n+1);
Pos(1,:)=[0,X];
Pos(2,:)=[1,X];
for i=1:n
    Pos(i+2,:)=[X(1:i),1-X(i),X(i+1:end)];
end