clear all
[n,v0,edges,r,x,ind_PVs] = TWONODES();
matpower_name = 'matpower_cases/case_two_nodes.m';

[n,v0,edges,r,x,ind_PVs] = SCE47_adapted(); % bus numbering: substation is node 1
matpower_name = 'matpower_cases/case_SCE47.m';

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
runs = 10^4; % 10^4 for thesis plots

range = [0:.001:.01]; % stepsize 0.5 for thesis plots
pointsP = zeros(1,length(range));
pointsS = zeros(1,length(range));
tic
for j=1:length(range)
    disp("progress bar "+ num2str(j)+ "/" +num2str(length(range)))
    mean = range(j);
    
    % twonode
%     mucon = mean*ones(1,n); 
%     mugen = 0*mean*ones(1,n); % 0 or .5
%     sigmacon = eye(n);
%     sigmagen = [1 1; 1 1];
%     sampleplus = max(mvnrnd(mucon,sigmacon,runs),0);
%     samplemin = max(mvnrnd(mugen,sigmagen,runs),0);
    
    % SCE47
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
end

probP = pointsP/runs;
probS = pointsS/runs;

figure
fontsize = 35;
linewidth = 3;
plot(range, probS,'LineWidth',linewidth)
hold on
plot(range, probP,'LineWidth',linewidth)
hold off
xlabel({'$m$'},'Interpreter','latex')
legend({'$\pi_S$','$\pi_P$'},'Interpreter','latex','FontSize',fontsize+5);
set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 20, 20], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
saveas(gca,'risk_sce47_1','epsc') %gcf

toc

%% THEORY (MODIFIED) (didnt work out in the end; simulations suffice)
% probmod = 1;
% for i = 1:2*n
%     probmod = probmod * normcdf((delta*v0-Amat(i,:)*mu')/norm(sqrtm(sigma)*Amat(i,:)',2)); % not the right formula
%     disp(normcdf((delta*v0-Amat(i,:)*mu')/norm(sqrtm(sigma)*Amat(i,:)',2)))
% end
% probmod

