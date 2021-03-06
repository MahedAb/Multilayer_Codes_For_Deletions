clear
%%parameters
Iteration=1000;
k=5; %Number of deletions
L_1=7; % Number of blocks
L_2=5; % Nyumber of chunk-strings 
ch=5;  % Chunck length
n_b=ch*L_2; %block length
Num_Lin_eq=40;  % Number of randon linear codes
N=ch*L_1*L_2;  %total code length
ee=0;  % Number of unresolved instances
%% Variables

Ave_Blocks=0;
Ave_Chunk=0;
Ave_Sec=0;
Ave_Lin=0;
Ave_List=0;
P_list=0;
P_e=0;
list_zero=0;

h = waitbar(0,'Please wait...');
%% Main loop of the Code
for i=1:Iteration
    Candidate_list=[];
    list_size=0;
    dis_flag=0;
    problem_flag=1;
    waitbar(i/Iteration)
    %Encoding
    [X,synf,syng,b,Y,R,r,A]=Code_Maker(ch,L_1,L_2,k,Num_Lin_eq);
    %First step of the decoding
    can=Tree_Search(X,ch,L_1,L_2,synf);
    % Second step of the decoding
    [Z,blocks]=Prim_fix(X,can,synf,ch*L_2);
    [u,~]=size(Z);
    Ave_Blocks=Ave_Blocks+u;
    
    %Steps 3 to 6
    for j=1:u
        
        %% If no deletion remains after step 2, go to step 6
        if sum(blocks(j,:))==0
            if ~Discard_constraints(Z(j,:),synf,syng,ch,L_1,L_2)
                
                Ave_Chunk=Ave_Chunk+1;
                Ave_Sec=Ave_Sec+1;
                if mod(A*Z(j,:)',2)==r
                    dis_flag=dis_flag+1;
                    Ave_Lin=Ave_Lin+1;
                    if list_size==0
                        Candidate_list=Z(j,:);
                        list_size=1;
                    elseif ~ismember(Z(j,:),Candidate_list,'rows')
                        list_size=list_size+1;
                        Candidate_list(list_size,:)=Z(j,:);
                    end
                end
            end
            
            %% Going to step 3 (when we have deletions)
        else
            [can1,s]=Chunk_Search_1(Z(j,:),ch,L_2,syng,blocks(j,:));
            Ave_Chunk=Ave_Chunk+s;
            for t=1:s
                [dis,X_1,H]=Second_Fix(Z(j,:),can1(:,:,t),blocks(j,:),synf,syng,ch,L_2);
                if dis==0
                    dis_flag=dis_flag+1;
                    Ave_Sec=Ave_Sec+1;
                    [X_2,S]=ready_Lin(X_1,H,ch);
                    [c,X_3]=Fourth_step(X_2,A,r);
                    % You can add a flag that show whther the answer was
                    % unique, which shows if you need to go for else
                    
                    %Step 4 of the decoding
                    if c==1;
                        Ave_Lin=Ave_Lin+1;
                        if ~Discard_constraints(X_3,synf,syng,ch,L_1,L_2)
                            % you can add substring criteria checking here
                            if list_size==0
                                Candidate_list=X_3;
                                list_size=1;
                            elseif ~ismember(X_3,Candidate_list,'rows')
                                list_size=list_size+1;
                                Candidate_list(list_size,:)=X_3;
                            end
                        else
                            % It means that possibly there are more than one solution for linear equations,
                            %and we got the wrong one, so here we guess one bit in every chunks with one or 2 deletions and then try to
                            %recover deletions using VT then solve the equations and see whether the
                            %result is consistent with the VTs.
                            flag=1;
                            elem=find(H);
                            for ww=elem'
                                if flag==0
                                    break
                                end
                                if H(ww)<3
                                    H1=H;
                                    col=ceil(ww/L_1);
                                    row=ww-(col-1)*L_1;
                                    % row and column of the chunk
                                    % that we are gussing bit in it
                                    ind=(row-1)*n_b+(col-1)*ch-sum(sum(H1(1:(row-1),:)))-sum(H1(row,1:(col-1)))+1;
                                    temp=X_1(ind:ind+ch-H1(ww)-1);
                                    H1(ww)=H1(ww)-1;
                                    %temp comprises the chunk
                                    pos=One_add(temp);
                                    % in next loop we add all
                                    % possible
                                    for rr=1:length(temp)+2
                                        X_2=[X_1(1:ind-1),pos(rr,:),X_1(ind+ch-H1(ww)-1:end)];
                                        [dis,X1_2,H2]=Second_Fix_1(X_2,H1,synf,syng,ch);
                                        if ~dis
                                            Ave_Sec=Ave_Sec+1;
                                            [X2_2,S]=ready_Lin(X1_2,H2,ch);
                                            [c,X_3]=Fourth_step(X2_2,A,r);
                                            if c
                                                Ave_Lin=Ave_Lin+1;
                                                if ~Discard_constraints(X_3,synf,syng,ch,L_1,L_2)
                                                    % you can add substring thing here
                                                    flag=0;
                                                    if list_size==0
                                                        Candidate_list=X_3;
                                                        list_size=1;
                                                        ee=ee+1;
                                                    elseif ~ismember(X_3,Candidate_list,'rows')
                                                        list_size=list_size+1;
                                                        Candidate_list(list_size,:)=X_3;
                                                        ee=ee+1;
                                                    end
                                                end
                                            end
                                        end
                                        if flag==0
                                            break
                                        end
                                    end
                                end
                                
                            end
                            
                            
                        end
                    end
                end
            end
        end
    end
    % If we list size is larger than one, we are making errors
    if list_size>1
        P_list=P_list+1;
        SOS=Candidate_list;
        P_e=P_e+sum(abs(Y-Candidate_list(1,:)))/N;
    % This means something is wrong! Check the deletion pattern.
    elseif list_size==0
        list_zero=list_zero+1;
        if problem_flag
            H1 
        end
        P_e=P_e+1;
    else
        P_e=P_e+sum(abs(Y-Candidate_list(1,:)))/N;
    end
    Ave_List=Ave_List+list_size;
end
close(h)

Ave_Blocks=Ave_Blocks/Iteration
Ave_Chunk=Ave_Chunk/Iteration
Ave_Sec=Ave_Sec/Iteration
Ave_Lin=Ave_Lin/Iteration

Ave_List=Ave_List/Iteration
P_list
list_zero
P_e=P_e/Iteration
Rate=1-(Num_Lin_eq/N)-log2(ch*L_2+1)/(ch*L_2)-log2(ch*L_1+1)/(ch*L_1)


