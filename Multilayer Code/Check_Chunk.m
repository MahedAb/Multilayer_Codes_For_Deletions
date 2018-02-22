function x=Check_Chunk(pointer,X,block,j,ch,l_2)
tt=0;
l_1=length(block);
n=ch*l_2;
s=zeros(1,l_1*ch);
for i=1:l_1
    B=(i-1)*n-sum(block(1:i-1));
    if block(i)==0
        s((i-1)*ch+1:i*ch)=X(B+(j-1)*ch+1:B+j*ch);
    else
        tt=tt+1;
        s((i-1)*ch+1:i*ch)=X(B+pointer(tt):B+pointer(tt)+ch-1);
    end
end
x=VT(s);
