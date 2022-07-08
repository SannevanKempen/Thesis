% initialize
% tic
clear all
tic
[n,v0,edges,r,x,ind_PVs] = TWONODES(); % bus numbering: substation is node 1
r = [0.01 0.01];
x = [0 0];
ind_PVs = [2 3];
savingfile = 'SOCP_vs_R_twonode.mat';

% [n,v0,edges,r,x,ind_PVs] = SCE47_adapted();
% x = zeros(1,n);
% savingfile = 'SOCP_vs_R_SCE47.mat';

% [n,v0,edges,r,x,ind_PVs] = TWENTYONENODES();
% ind_PVs = 2:22;
% savingfile = 'SOCP_vs_R_22.mat';
% [n,v0,edges,r,x,ind_PVs] = SCE56();
% savingfile = 'SOCP_vs_R_SCE47.mat';

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


%% OPF (SOCP + MOD)
n = n-1; % bus numbering: substation is node 0
pconub = 10*ones(1,n); % 10 for twonode; .01 for SCE47
qconub = tan(acos(1))*pconub;
pgenub = zeros(1,n);
qgenub = zeros(1,n);

range = [0:2:2]; % [0:0.1:15] for twonode; [0:.005:.15] for SCE47
optval_mod = zeros(1,length(range));
optval_socp = zeros(1,length(range));
powers = zeros(length(range),12*n);
for i = 1:length(range)
    ind = range(i)
    pgenub(ind_PVs - 1) = ind; % ind for twonode; 
    qgenub = tan(acos(1))*pgenub;
    [pcon_mod,qcon_mod,pgen_mod,qgen_mod,pcon_socp,qcon_socp,pgen_socp,qgen_socp,volt_socp,W,WY,time_mod,time_socp] = OPF_SOCP_and_R(n,Ymat,Amat,v0,delta,pconub,qconub,pgenub,qgenub);
    powers(i,:) = [pconub qconub pgenub qgenub pcon_mod qcon_mod pgen_mod qgen_mod pcon_socp(2:end) qcon_socp(2:end) pgen_socp(2:end) qgen_socp(2:end)];
    optval_mod(i) = sum(log(pcon_mod));
    optval_socp(i) = sum(log(pcon_socp(2:end)));
end
toc

%% SAVE
save(savingfile,'optval_mod','optval_socp','range','pconub');

%% ANALYSIS 
% error1 = (powers(:,8*n+1) - powers(:,4*n+1))./powers(:,8*n+1)
% error2 = (powers(:,8*n+2) - powers(:,4*n+2))./powers(:,8*n+2)

%% PLOT

% load('SOCP_vs_R_twonode.mat');
% load('SOCP_vs_R_SCE47.mat');
figure
fontsize = 45;
linewidth = 3;
plot(range, optval_socp,range, optval_mod,'LineWidth',linewidth)
% plot(range,(optval_socp-optval_mod)./-optval_socp,'LineWidth',linewidth)
grid on
% ylim([0 0.02]);
% xlim([0 range(end)]);

xlabel({'$g$'},'Interpreter','latex','FontSize',fontsize+5)
legend({'SDP','P'},'Location','northwest','FontSize',fontsize+5);
set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 20, 20], 'PaperUnits', 'Inches', 'PaperSize', [30, 30])
saveas(gca,'Accuracy','epsc') %gcf
% saveas(gca,'AccuracyRelError','epsc') %gcf


%% SAVING RESULTS 
% save(savingfile,'[pcon_mod pcon_socp(2:end)]')
