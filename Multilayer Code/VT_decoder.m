function [X,b,p]=VT_decoder(x,a)
n=length(x)+1;
y=0;
w=sum(x);
for i=1:n-1
    y=y+i*x(i);
end

y=mod(y,n+1);
D=mod(a-y,n+1);

if D<=w
    b=0;
    R_1=D;
    i=n-1;
    while sum(x(i:end))<R_1
        i=i-1;
    end
    i=i+(D==0);
    p=i;
    X=[x(1:i-1),0,x(i:end)];
end
if D>w
    b=1;
    L_0=D-w-1;
    y=mod(x+1,2);
    i=1;
    while sum(y(1:i))<L_0
        i=i+1;
    end
    p=i+1-(L_0==0);
    i=i-(L_0==0);
    X=[x(1:i),1,x(i+1:end)];
end
