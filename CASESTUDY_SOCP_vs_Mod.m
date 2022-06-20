% initialize
clear all
% [n,v0,edges,r,x,ind_PVs] = TWONODES(); % bus numbering: substation is node 1
[n,v0,edges,r,x,ind_PVs] = SCE47_adapted();

delta = 0.1;

%% CONSTRUCTION Y/Z/A MATRIX
Ymat = zeros(n,n);
for i=1:length(edges)
    Ymat(edges(i,1),edges(i,2)) = -1/(r(i));
    Ymat(edges(i,2),edges(i,1)) = -1/(r(i));
end
for i=1:n
    Ymat(i,i) = -sum(Ymat(i,:));
end
Zmat = inv(Ymat(2:end,2:end));
Amat = [1/((1-delta)*v0)*Zmat -1/((1+delta)*v0)*Zmat; -1/((1+delta)*v0)*Zmat 1/((1-delta)*v0)*Zmat];

%% BOUNDS
n = n-1; % bus numbering: substation is node 0
pconub = 100*ones(1,n);
pgenub = zeros(1,n);
pgenub(ind_PVs-1) = .1;

%% SOCP
[pcon_SOCP,pgen_SOCP,W] = OPF_SOCP(n,Ymat,v0,delta,pconub,pgenub);

%% MODIFIED
[pcon_mod,pgen_mod] = OPF_modified(n,Amat,v0,delta,pconub,pgenub);

%% ANALYSIS
optval_SOCP = sum(arrayfun(@(x) log(x),pcon_SOCP(2:end)))
opt_val_mod = sum(arrayfun(@(x) log(x), pcon_mod))
pcon_SOCP(2:end) - pcon_mod