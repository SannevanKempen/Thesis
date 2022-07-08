clear all
lambda = .5;
mu = 1;

t = [.1 : .1 : .9]; %[0.01:0.01:1];
nrRuns = 50;
EZ1 = zeros(nrRuns,length(t));
EZ2 = zeros(nrRuns,length(t));
Epower1 = zeros(nrRuns,length(t));
Epower2 = zeros(nrRuns,length(t));
for i = 1:length(t)
    delta = t(i)
    tic
    for j=1:nrRuns
        [EZ1(j,i),EZ2(j,i),Epower1(j,i),Epower2(j,i),runningEZ1,runningEZ2] = Simulation(delta,lambda,mu); 
    end
    toc
end
save('Simulationresults.mat','EZ1','EZ2','Epower1','Epower2')


EZ1SEM = std(EZ1)/sqrt(nrRuns); % for confidence intervals
EZ2SEM = std(EZ2)/sqrt(nrRuns);
EZ1Mean = mean(EZ1);
EZ2Mean = mean(EZ2);
CI95 = tinv([0.025 0.975], nrRuns-1);  
EZ1CI95 = bsxfun(@times, EZ1SEM, CI95(:));
EZ2CI95 = bsxfun(@times, EZ2SEM, CI95(:));

r1 = .01;
r2 = .015;
s = [.001 : .001 : 1];
rhoLDF1 = 2./(1-(1-s).^2)*lambda*r1*mu;
rhoR1 = 1./((1-s).*s)*lambda*r1*mu;

rhoLDF2 = 2./(1-(1-s).^2)*lambda*r2*mu;
rhoR2 = 1./((1-s).*s)*lambda*r2*mu;


% figure
% fontsize = 45;
% linewidth = 3;
% plot(s,rhoR1,s,rhoLDF1,'LineWidth',linewidth)
% xlim([0.05 .95])
% % ylim([0 1])
% xlabel({'$\Delta$'},'Interpreter','latex','FontSize',fontsize+5)
% ylabel({'$\rho_1$'},'Interpreter','latex','FontSize',fontsize+5)
% legend({'$\textrm{P}$', '$\textrm{LDF}$','\textrm{DF}'},'Interpreter','latex','FontSize',fontsize+5,'Location','southwest');
% set(gca,'FontSize',fontsize)
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 20, 20], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
% saveas(gca,'rho1','epsc') %gcf
% 
% figure
% fontsize = 45;
% linewidth = 3;
% plot(s,rhoR2,s,rhoLDF2,'LineWidth',linewidth)
% xlim([0.05 .95])
% % ylim([0 1])
% xlabel({'$\Delta$'},'Interpreter','latex','FontSize',fontsize+5)
% ylabel({'$\rho_2$'},'Interpreter','latex','FontSize',fontsize+5)
% legend({'$\textrm{P}$', '$\textrm{LDF}$','\textrm{DF}'},'Interpreter','latex','FontSize',fontsize+5,'Location','southwest');
% set(gca,'FontSize',fontsize)
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 20, 20], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
% saveas(gca,'rho2','epsc') %gcf


EZLDF = (rhoLDF1 + rhoLDF2)./(ones(1,length(rhoLDF1))-(rhoLDF1 + rhoLDF2)); % mean nr of EVs
EZR = (rhoR1 + rhoR2)./(ones(1,length(rhoR1))-(rhoR1 + rhoR2));
EZDF = EZ1Mean+EZ2Mean;

% 
figure
fontsize = 45;
linewidth = 3;
plot(s,EZR,s,EZLDF,'LineWidth',linewidth)
hold on 
plot(t, EZDF,'x','MarkerSize',25,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
% hold on 
% errorbar(t,EZDF,EZ1CI95(2,:),'o','MarkerSize',10)
hold off
grid on;
xlim([0.05 .95])
% ylim([0 1])
xlabel({'$\Delta$'},'Interpreter','latex','FontSize',fontsize+5)
ylabel({'EVs'},'Interpreter','latex','FontSize',fontsize+5)
legend({'$\textrm{P}$', '$\textrm{LDF}$','\textrm{DF}'},'Interpreter','latex','FontSize',fontsize+5,'Location','southwest');
set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 20, 20], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
saveas(gca,'ProdQueueEZ','epsc') %gcf
