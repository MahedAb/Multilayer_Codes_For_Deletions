%This function computes VT syndrome of x
function y=VT(x)
y=0;
n=length(x);
for i=1:n
    y=y+i*x(i);
end
y=mod(y,n+1);

% n=length(x);
% i=(1:n);
% y=mod(sum(i.*x),n+1);