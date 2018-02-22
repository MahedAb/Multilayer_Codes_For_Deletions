function [Discard,X,H]=Second_Fix_1(X,H,synf,syng,t)
Discard=0;
[L_1,L_2]=size(H);
n=t*L_2;
%Discard=1 shows that the candidate get discarded
blo=sum(H,2)';

s=zeros(1,n-1);
while (sum(sum(H,2)==1)>0 || sum(sum(H)==1)>0)&& Discard==0;
    %First for, for g functions with single deletions
    for i=find(sum(H)==1)
        point=1;
        j=find(H(:,i));
        for r=1:L_1
            pointb=(r-1)*n+(i-1)*t-sum(H(r,1:i-1))-sum(blo(1:r-1))+1;
            if r==j
                s(point:point+t-2)=X(pointb:pointb+t-2);
                point=point+t-1;
            else
                
                s(point:point+t-1)=X(pointb:pointb+t-1);
                point=point+t;
            end
        end

        [~,b,kk]=VT_decoder(s,syng(i));
   %     kk

        if (kk>j*t)
            while s(kk-1)==b && kk>j*t
                kk=kk-1;
            end
            if (kk>j*t)
                Discard=1;
                break;
            end
        elseif kk<(j-1)*t+1
            while s(kk)==b && kk<(j-1)*t+1
                kk=kk+1;
            end
            if kk<(j-1)*t+1
                Discard=1;
                break;
            end
        end
  %      kk
        kk=kk-(j-1)*t;
        pos=n*(j-1)+(i-1)*t-sum(H(j,1:i-1))-sum(blo(1:j-1))+kk;
        X=[X(1:pos-1),b,X(pos:end)];
        H(j,i)=0;
        blo(j)=blo(j)-1;
    end
 %   Discard
    if Discard==1 || sum(sum(H))==0
        break;
    end
    % Second for, for blocks that now have single deletion
    tt=n*L_1-sum(blo);
    for i=find((sum(H,2)==1)')
%         if sum(sum(H,2)==1)==0
%             break;
%         end
        
        
        blo(i)=0;
        j=find(H(i,:));
        H(i,:)=0;
        
        pointer=(i-1)*n-sum(blo(1:i-1));
        
        [S,b,kk]=VT_decoder(X(pointer+1:pointer+n-1),synf(i));
        if (kk>j*t)
            %Are you sure about the pointer+kk vs pointer+kk+1 
            while X(pointer+kk-1)==b && kk>j*t
                kk=kk-1;
            end
            if (kk>j*t)
                Discard=1;
                break;
            end
        elseif kk<(j-1)*t+1
            while X(pointer+kk)==b && kk<(j-1)*t+1
                kk=kk+1;
            end
            if kk<(j-1)*t+1
                Discard=1;
                break;
            end
        end
        X(1:tt+1)=[X(1:pointer),S,X(pointer+n:tt)];
        tt=tt+1;
    end
end
%X=X(1:n*L_1);