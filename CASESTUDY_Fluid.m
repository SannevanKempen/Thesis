% initialize

tic
clear all
% [n,v0,edges,r,x,ind_PVs] = TWONODES(); % bus numbering: substation is node 1
[n,v0,edges,r,x,ind_PVs] = SCE47_adapted();
savingfile = 'ChargingProbSCE47.mat';
% [n,v0,edges,r,x,ind_PVs] = SCE56();
% savingfile = 'ChargingProbSCE56.mat';

delta = 0.1;

%% CONSTRUCTION Y/Z/A MATRIX
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
% Xmat = zeros(n-1,n-1);
% Amat = [1/((1-delta)*v0)*Zmat -1/((1+delta)*v0)*Zmat; -1/((1+delta)*v0)*Zmat 1/((1-delta)*v0)*Zmat];
Amat= [1/((1-delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat   -1/((1+delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat;
    1/((1-delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat   -1/((1+delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat;
    -1/((1+delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat    1/((1-delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat;
    -1/((1+delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat    1/((1-delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat];

%% PARAMETERS
n = n-1; % bus numbering: substation is node 0
pconub = 10^3*ones(1,n); % practically unbounded
pgenub = zeros(1,n);
pgenub(ind_PVs - 1) = 0;
qgenub = tan(acos(1))*pgenub;
lambda = 1/46*ones(1,n); %1/n*ones(1,n);
mu = 1*ones(1,n);
nu = .5*ones(1,n);
K = 1*ones(1,n);
c = 1;

%% OPF (SOCP + MOD)
[pcon_mod,pgen_mod,qgen_mod,pcon_socp,pgen_socp,qgen_socp,volt_fluid] = OPF_fluid(n,Ymat,Amat,v0,delta,pconub,pgenub,qgenub,lambda,mu,nu,K,c);
power = [pcon_socp(2:end); pcon_mod];
%% ANALYSIS 
lambdaK = Erlangloss(n,K,lambda,nu);
z = zeros(2,n);
probs = zeros(2,n);
for type = 1:2
for i = 1:n
    if lambda(i)>0
        z(type,i) = (lambdaK(i)-mu(i) * power(type,i))/nu(i);
        probs(type,i) = 1 - z(type,i)*nu(i)/lambdaK(i);
    end
end
end
power
z

xi = sum(power(1,:))
lam = sum(z(1,:))
omega = 1-sum(z(1,:))/sum(lambdaK./nu)

xi = sum(power(2,:));
lam = sum(z(2,:));
omega = 1-sum(z(2,:))/sum(lambdaK./nu);

probs;

sum(arrayfun(@(x) log(x),power(1,:)));
sum(arrayfun(@(x) log(x),power(2,:)));


%% SAVING RESULTS FOR COLOR PLOTS (WITH PYTHON SCRIPT IN FOLDER 'TikZ_Color_SCE')
p = probs(1,:);
save(savingfile,'p')

toc