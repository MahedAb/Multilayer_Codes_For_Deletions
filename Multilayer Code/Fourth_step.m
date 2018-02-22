function [c,X]=Fourth_step(X,A,r)

c=1; % if c=0 then some of constraints are not satisfied so candidate should be discarded
flag=1; % Internal use
%flag_1=1; % if falg_1=0 then we failed to solve the equations
[m,~]=size(A);
I=1:m;

%checking other constraints
while c==1 && flag
    flag=0;
    for i=I
        if sum(X(A(i,:))<0)==0
            
            if mod(sum(X(A(i,:))),2)~=r(i)
                c=0;
            else
                I=I(I~=i);
            end
        elseif sum(X(A(i,:))<0)==1
            
            flag=1;
            j=find(X.*A(i,:)<0);
            X(j)=mod(sum(X(A(i,:)))+r(i)+1,2);
        end
    end
end


%Solving the linear equations
if c==1 && sum(X<0)>0
    A1=A(:,X<0);
    A2=A(:,X>=0);
    r2=mod(A2*X(X>=0)',2);
    r1=mod(r-r2,2);
    [x,flag_1]=gflineq(A1,r1);
    if flag_1==1
        X(X<0)=x;
        s=mod(A*X',2);
        if s~=r
            c=0;
        end
    else
        c=0;
    end
end
