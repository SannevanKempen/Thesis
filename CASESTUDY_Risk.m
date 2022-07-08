clear all

% [n,v0,edges,r,x,ind_PVs] = TWONODES();
% r = [0.01 0.01];
% x = [0 0];
% ind_PVs = [];
% matpower_name = 'matpower_cases/case_two_nodes.m';
% % savingfile = 'Risk_twonode_g0.mat';
% savingfile = 'Risk_twonode_g0.5.mat';

[n,v0,edges,r,x,ind_PVs] = SCE47_adapted(); % bus numbering: substation is node 1
x = zeros(1,46);
matpower_name = 'matpower_cases/case_SCE47.m';
% savingfile = 'Risk_SCE47_g0.mat';
savingfile = 'Risk_SCE47_g3.mat';

% [n,v0,edges,r,x,ind_PVs] = TWENTYONENODES();
% ind_PVs = [14 22];
% matpower_name = 'matpower_cases/case_21.m';

delta = 0.1;

%% CONSTRUCTION Z/A MATRIX
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

mpc = loadcase(matpower_name);
mpc.gen(1, 6) = v0; % set correct values for v0,r and x in mpc
for j=1:n-1
    mpc.branch(j, 3) = r(j);
    mpc.branch(j, 4) = x(j);
end
%% SIMULATIONS (MODIFIED + ORIGINAL OPF)
n = n-1;
range = logspace(-3,-2,30); % 30 points. Interval [0,1] for twonode g=0; [0,1.2] g=0.5; [-3,-2] SCE-47 g=0; [-3,-2] g=3
pointsP = zeros(1,length(range));
pointsS = zeros(1,length(range));

for j=1:length(range)
    tic
    disp("progress bar "+ num2str(j)+ "/" +num2str(length(range)))
    mean = range(j);
    if j < length(range)/4
        runs = 5*10^4; %5*10^4;
    else
        runs = 5*10^4; %5*10^3;
    end
    % twonode
%     mucon = mean*ones(1,n); 
%     mugen = 0.5*mean*ones(1,n); % 0 or .5
%     sigmacon = eye(n);
%     sigmagen = ones(n,n);
%     sampleplus = max(mvnrnd(mucon,sigmacon,runs),0);
%     samplemin = max(mvnrnd(mugen,sigmagen,runs),0);
    
%     % SCE47
    mucon = mean*ones(1,n);
    sigmacon = 0.0001*eye(n);
    sampleplus = max(mvnrnd(mucon,sigmacon,runs),0);
    nPVs = length(ind_PVs);
    samplemin = zeros(runs,n);
    mugen = 3*mean*ones(1,nPVs); % 0 or 3
    sigmagen = 0.0001*ones(nPVs,nPVs);
    samplemin(:,ind_PVs - 1) = max(mvnrnd(mugen,sigmagen,runs),0);
    
    for i=1:runs
        if Amat*[sampleplus(i,:) samplemin(i,:)]' <= delta*v0
            pointsP(j) = pointsP(j) + 1;
        end
        
        for k=1:n
            mpc.bus(k+1, 3) = sampleplus(i,k) - samplemin(i,k);
        end
        [results,succes] = runpf(mpc, mpoption('verbose',0,'out.all',0));
        if succes == 1
            nr_of_violations = sum(results.bus(:,8) < (1-delta)*v0) + sum(results.bus(:,8) > (1+delta)*v0);
            if nr_of_violations == 0
                pointsS(j) = pointsS(j) + 1;
            end
        end
    end
    toc
end

probP = pointsP/runs;
probS = pointsS/runs;
save(savingfile,'probP','probS','range','runs');

clear all
load('Risk_twonode_g0.mat');
% load('Risk_twonode_g0.5.mat');
% load('Risk_SCE47_g0.mat');
% load('Risk_SCE47_g3.mat');
probS = 1-probS;
probP = 1-probP;

figure
fontsize = 30;
linewidth = 2;
loglog(range, probP,range, probS,'LineWidth',linewidth)
grid on
ylim([.01 1])
yticks([0.05 0.1 0.5 1])
% yticks([.4:.1:1])
% xlim([1 10])
xlabel({'$m$'},'Interpreter','latex','FontSize',fontsize+5)
legend({'$\pi_P$','$\pi_S$'},'Location','northwest','Interpreter','latex','FontSize',fontsize+5);
set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12, 20], 'PaperUnits', 'Inches', 'PaperSize', [10, 20])
saveas(gca,'Risk_twonode_g0','epsc') %gcf



%% THEORY (MODIFIED) (didnt work out in the end; simulations suffice)
% probmod = 1;
% for i = 1:2*n
%     probmod = probmod * normcdf((delta*v0-Amat(i,:)*mu')/norm(sqrtm(sigma)*Amat(i,:)',2)); % not the right formula
%     disp(normcdf((delta*v0-Amat(i,:)*mu')/norm(sqrtm(sigma)*Amat(i,:)',2)))
% end
% probmod

