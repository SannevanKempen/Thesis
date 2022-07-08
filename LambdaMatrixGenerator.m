% Script to determine power allocation function in network of the form 0 -- 1 -- 2 with possible RES on node 2
% Used as input for mathematica script ('MarkovChain.m') to compute stat distr of the 4D Markov chain
% Optimization properties
clear;
tic
range = [.1 : .1 : .9];
for m=1:length(range)
    vars = {'delta','Ymat','Zmat','Rmat','Xmat','Amat','LambdaMatrix','temp','n','alc_mod','alc_socp','alc_ldf','socp'};
    clear(vars{:})
    n = 3;
    edges = [1 2; 2 3];
    r = [0.01 0.005];
    x = [0 0];
    u = [1 1];
    v0 = 1;
    K = 20;
    P2gen=0;
    Q2gen=0;
    
    delta = range(m)
    Ymat = zeros(n,n);
    for i=1:length(edges)
        Ymat(edges(i,1),edges(i,2)) = -1/(r(i)+1i*x(i));
        Ymat(edges(i,2),edges(i,1)) = -1/(r(i)+1i*x(i));
    end
    for i=1:n
        Ymat(i,i) = -sum(Ymat(i,:));
    end
    Zmat = inv(Ymat(2:end,2:end));
    Rmat = real(Zmat);
    Xmat = imag(Zmat);
    Amat= [1/((1-delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat   -1/((1+delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat;
        1/((1-delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat   -1/((1+delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat;
        -1/((1+delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat    1/((1-delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat;
        -1/((1+delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat    1/((1-delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat];
    
    
    
    LambdaMatrix = zeros(K+2,K+2,9); % need larger matrix first/last row/column will not be used in mathematica calculations anyway
    temp = 1;
    n = 2;
    tic
    for z1=0:K+1
        for z2=0:K+1
            if z1+z2>0 % ignore case z1=z2=0
                z = [z1 z2];
                [alc_mod,alc_socp,alc_ldf] = OPF_SOCP_and_R_and_LDF(n,Ymat,Zmat,Amat,v0,delta,u,P2gen,Q2gen,z);
                LambdaMatrix(z1+1,z2+1,:) = [0 alc_mod 0 alc_ldf alc_socp];
                %             % for fluid:
                % LambdaMatrix(z1+1,z2+1,:) = OPT(I,A,R,X,u,VU,VU,AC_bool,P2gen,Q2gen,z);
                
                %             % for modified:
                %             [pmin,pplus,dual1,dual2] = OPF_Casestudy_twonodes(n,Amat,u(2:end),delta,v0,z(2:end),[0 P2gen]);
                %             p_modified = pplus-pmin;
                %             LambdaMatrix(z1+1,z2+1,:) = [0 p_modified];
                
                disp(['progress bar ', num2str(temp),'/', num2str((K+2)^2)])
                temp = temp + 1;
            end
        end
    end
    toc
    LambdaMatrix;
    
    mod = LambdaMatrix(:,:,1:3);
    ldf = LambdaMatrix(:,:,4:6);
    socp = LambdaMatrix(:,:,7:9);
    % save(['LambdaMatrices/mod' num2str(K) '_' num2str(P2gen) '.mat'],'mod');
    % save(['LambdaMatrices/ldf' num2str(K) '_' num2str(P2gen) '.mat'],'ldf');
    % save(['LambdaMatrices/socp' num2str(K) '_' num2str(P2gen) '.mat'],'socp');
    
%         save(['LambdaMatrices/LambdaMatrixSimulation' num2str(delta) '.mat'],'socp');
    % test=load('LambdaMatrices/test0.2.mat')
end