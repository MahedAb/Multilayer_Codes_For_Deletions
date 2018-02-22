%This function gives us possible candidates for block boundaries and inputs
%are X received codeword and t ,k parameters of code and also syndroms
function [Can]=Tree_Search(X,t,L_1,L_2,fsyn)

n=t*L_2;
k=t*L_1*L_2-length(X);
Max=nchoosek(L_1+k-1,k);
Del=k*ones(Max,1);
pointer=ones(Max,1);
dead=zeros(Max,1);
H=zeros(Max,L_1);
Num=1;

for i=1:L_1-1
    branch=Num;
    for j=1:branch
        if ~dead(j)
            if VT(X(pointer(j):pointer(j)+n-1))==fsyn(i)
                H(j,i)=0;
                
                H(Num+1:Num+1+Del(j)-1,1:i-1)=repmat(H(j,1:i-1),Del(j),1);
                for tt=2:Del(j)
                    Del(Num+tt-1)=Del(j)-tt;
                    pointer(Num+tt-1)=pointer(j)+n-tt;
                    H(Num+tt-1,i)=tt;
                end
                Num=Num+Del(j)-1+(Del(j)==0);
                pointer(j)=pointer(j)+n;
            else
                if Del(j)<1
                    dead(j)=1;
                else
                    H(j,i)=1;
                    H(Num+1:Num+1+Del(j)-1,1:i-1)=repmat(H(j,1:i-1),Del(j),1);
                    for tt=2:Del(j)
                        Del(Num+tt-1)=Del(j)-tt;
                        pointer(Num+tt-1)=pointer(j)+n-tt;
                        H(Num+tt-1,i)=tt;
                    end
                    Num=Num+Del(j)-1;
                    Del(j)=Del(j)-1;
                    pointer(j)=pointer(j)+n-1;
                end
            end
        end
        
    end
end
for j=1:Num
    if ~dead(j)
        if Del(j)==0
            if VT(X(pointer(j):pointer(j)+n-1))~=fsyn(L_1)
                dead(j)=1;               
            end            
        elseif Del(j)~=n-(length(X)-pointer(j)+1)
            c=j
            12334
        else
            H(j,L_1)=Del(j);
            Del(j)=0;
        end
    end
end
H=H(1:Num,:);
dead=dead(1:Num);
dead=dead<0.5;
Can=H(dead,:);