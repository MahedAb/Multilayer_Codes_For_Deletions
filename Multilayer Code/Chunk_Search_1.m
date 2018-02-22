function [Can,survive]=Chunk_Search_1(X,t,L_2,syn,block)

s=nnz(block);
L_1=length(block);
M=(prod(block+1))^2;
Del=ones(s,M);
Del(:,1)=nonzeros(block);

%Since input string have some additional zeros at the end
X=X(1:t*L_1*L_2-sum(block));

pointer=ones(s,M);
dead=zeros(M,1);
H=zeros(s,L_2,M);
Num=1;

if sum(block)>0
    for i=1:L_2-1
        
        u=sum(dead==1);
        H=H(:,:,dead<1);
        pointer=pointer(:,dead<1);
        
        Del=Del(:,dead<1);
        dead=dead(dead<1);
        
        Num=Num-u;
        branch=Num;
        for j=1:branch
            if ~dead(j)
                if Check_Chunk(pointer(:,j),X,block,i,t,L_2)==syn(i)
                    
                    P=poss_matrix(Del(:,j),1,t);
                    H(:,i,j)=P(1,:)';
                    [len,~]=size(P);
                    
                    H(:,1:i-1,Num+1:Num+1+len-1)=repmat(H(:,1:i-1,j),1,1,len);
                    for tt=2:len
                        Del(:,Num+tt-1)=Del(:,j)-P(tt,:)';
                        pointer(:,Num+tt-1)=pointer(:,j)+t-P(tt,:)';
                        H(:,i,Num+tt-1)=P(tt,:)';
                        dead(Num+tt-1)=0;
                    end
                    
                    Num=Num+len-1;
                    pointer(:,j)=pointer(:,j)+t;
                    
                else
                    if sum(Del(:,j))<1
                        dead(j)=1;
                    else
                        
                        
                        P=poss_matrix(Del(:,j),0,t);
                        H(:,i,j)=P(1,:)';
                        [len,~]=size(P);
                        
                        H(:,1:i-1,Num+1:Num+1+len-1)=repmat(H(:,1:i-1,j),1,1,len);
                        for tt=2:len
                            Del(:,Num+tt-1)=Del(:,j)-P(tt,:)';
                            pointer(:,Num+tt-1)=pointer(:,j)+t-P(tt,:)';
                            H(:,i,Num+tt-1)=P(tt,:)';
                            dead(Num+tt-1)=0;
                        end
                        Num=Num+len-1;
                        pointer(:,j)=pointer(:,j)+t-P(1,:)';
                        Del(:,j)=Del(:,j)-P(1,:)';
                        
                    end
                end
            end
            
        end
%         Num
%         len
%         branch
%         length(pointer)
% 
%         length(dead)
        %     i
        %     dead(1:Num)
        %     H(1,:,1:Num)
    end
    for j=1:Num
        if ~dead(j)
            if sum(Del(:,j))==0
                if  Check_Chunk(pointer(:,j),X,block,L_2,t,L_2)~=syn(L_2)
                    dead(j)=1;
                end
                %         elseif Del(j)~=n-(length(X)-pointer(j)+1)
                %             c=j
            else
                H(:,L_2,j)=Del(:,j);
                Del(:,j)=0;
            end
        end
    end
end
H=H(:,:,1:Num);
dead=dead(1:Num);
dead=dead<0.5;
survive=sum(dead)-(sum(block)==0);
Can=H(:,:,dead);