%This function gives the number of deletions in each block based on blocks
%that had deletions

function x=indic_b(b,L)
x=zeros(1,L);
for i=1:length(b)
    x(b(i))=x(b(i))+1;
end