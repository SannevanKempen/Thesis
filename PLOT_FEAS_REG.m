clear all
tic
warning off;
warning('off','MATLAB:singularMatrix')

% TWO NODES
% [n,v0,edges,r,x,ind_PVs] = TWONODES();
n = 3;
v0 = 1;
edges = [1 2; 2 3];
delta = 0.1;
r = [0.01 0.01];
% x = [0.01 0.01];
x = [0.001 0.001];
x = 0.01*tan(acos(.9))*ones(1,2);
matpower_name = 'matpower_cases/case_two_nodes.m';

%% matpower
mpc = loadcase(matpower_name);
mpc.gen(1, 6) = v0; % set correct values for v0,r and x in mpc
for j=1:n-1
    mpc.branch(j, 3) = r(j);
    mpc.branch(j, 4) = x(j);
end

%% Y/Z/R/X/A mat
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
% Amat = [1/((1-delta)*v0)*Rmat 1/((1-delta)*v0)*Xmat -1/((1+delta)*v0)*Rmat -1/((1+delta)*v0)*Xmat;
%     -1/((1+delta)*v0)*Rmat -1/((1+delta)*v0)*Xmat 1/((1-delta)*v0)*Rmat 1/((1-delta)*v0)*Xmat;
%     1/((1-delta)*v0)*Xmat -1/((1+delta)*v0)*Rmat -1/((1+delta)*v0)*Xmat 1/((1-delta)*v0)*Rmat;
%     -1/((1+delta)*v0)*Xmat 1/((1-delta)*v0)*Rmat 1/((1-delta)*v0)*Xmat -1/((1+delta)*v0)*Rmat]
Amat= [1/((1-delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat   -1/((1+delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat;
    1/((1-delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat   -1/((1+delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat;
    -1/((1+delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat    1/((1-delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat;
    -1/((1+delta)*v0)*Rmat-1/((1+delta)*v0)*Xmat   -1/((1+delta)*v0)*Xmat+1/((1-delta)*v0)*Rmat    1/((1-delta)*v0)*Rmat+1/((1-delta)*v0)*Xmat    1/((1-delta)*v0)*Xmat-1/((1+delta)*v0)*Rmat];

%% FEASIBILITY PLOTS

% savingfile = 'PolyPlotx0.001p0r0.01d0.1zoom.mat';
% savingfile = 'PolyPlotx0.001pc3r0.01d0.1zoom.mat';
% savingfile = 'PolyPlotx0.001pg3r0.01d0.1zoom.mat';

pmin = -50;
qmin = -50;
pmax = 10;
qmax = 40;

step = 1; % for thesis: .05 for large plots and 0.01 for zoom (takes loooong time)
prange = [pmin:step:pmax];
qrange = [qmin:step:qmax];

S_feas = zeros(length(prange),length(qrange));
R_feas = zeros(length(prange),length(qrange));
B_feas = zeros(length(prange),length(qrange));

for i=1:length(prange)
    disp(prange(i)) % manual progress bar
    for j=1:length(qrange)
        pcon = [0; 0]; % 3
        pgen = [0; 0]; % 3
        qcon = [0; 0]; % tan(acos(.9))*pcon(1)
        qgen = [0; 0]; % tan(acos(.9))*pgen(1)
        if prange(i) >= 0
                        pcon(2) = prange(i);
%             pcon(1) = prange(i);
        else
                        pgen(2) = -prange(i);
%             pgen(1) = -prange(i);
        end
        if qrange(j) >= 0
                        qcon(2) = qrange(j);
%             pcon(2) = qrange(j);
        else
                        qgen(2) = -qrange(j);
%             pgen(2) = -qrange(j);
        end
        p = pcon-pgen;
        q = qcon-qgen;
        s = p+1i*q;
        
        if Amat*[pcon; qcon; pgen; qgen] <= delta*v0
            R_feas(i,j) = 1;
        end
        
        
        if norm(Ymat(2:end,2:end)\diag(conj(s)),inf) <= (delta-delta^2)*v0^2
            B_feas(i,j) = 1;
        end
        for k=1:n-1
            mpc.bus(k+1, 3) = p(k);
            mpc.bus(k+1, 4) = q(k);
        end
        [results,succes] = runpf(mpc,mpoption('verbose',0,'out.all',0));
        if succes == 1
            nr_of_violations = sum(results.bus(:,8) < (1-delta)*v0) + sum(results.bus(:,8) > (1+delta)*v0);
            if nr_of_violations == 0
                S_feas(i,j) = 1;
            end
        end
    end
end
sum(sum(S_feas-R_feas<0))
toc

% 'PolyPlotx0q0r0.01d0.1.mat'
% 'PolyPlotx0q0r0.01d0.3.mat'
% 'PolyPlotx0q0r0.05d0.1.mat'

% 'PolyPlotx0.001p0r0.01d0.1.mat'
% 'PolyPlotx0.001pc3r0.01d0.1.mat'
% 'PolyPlotx0.001pg3r0.01d0.1.mat'
filename = 'PolyPlotx0.001p0r0.01d0.1.mat';
load(filename)
% plotting
result = S_feas + 1*R_feas + 3*B_feas;%
fontsize = 35;
mymap1 = [1 1 1;1 0 0;0 0 1;1 1 0; 1 1 0;0 1 0];
h = figure;
h = pcolor(prange,qrange,result');
colormap(mymap1)
set(gca, 'CLim', [0, 5]) 
set(h, 'EdgeColor', 'none');
% axes and grid
%     xlim([-5 5]);
%     ylim([-10 10]);
xlim([pmin pmax]);
ylim([qmin qmax]);
xl = xlim;
yl = ylim;
hold on; % add grid manually (could not fix this in another way unfortunately)
for row = yl(1) : 5 : yl(end) % p2range(1) : 1 : p2range(end)
    line([xl(1), xl(end)], [row, row], 'Color', 'black'); % line([p1range(1), p1range(end)], [row, row], 'Color', 'black');
end
for col = xl(1) : 5 : xl(end) % p1range(1) : 1 : p1range(end)
    line([col, col], [yl(1), yl(end)], 'Color', 'black'); % line([col, col], [p2range(1), p2range(end)], 'Color', 'black');
end
hold off
lbs = {'$S$','$P$','$B$'};
cols = [patch(NaN,NaN,mymap1(2,:)) patch(NaN,NaN,mymap1(3,:)) patch(NaN,NaN,mymap1(4,:))];
legend(cols, lbs,'Interpreter','latex','FontSize',fontsize,'Location','northeast');
%     xlabel('$p_1$','Interpreter','latex')
%     ylabel('$p_2$','Interpreter','latex')
xlabel('$p_2$','Interpreter','latex')
ylabel('$q_2$','Interpreter','latex')
set(gca,'FontSize',fontsize)
xl = get(gca,'XLabel'); % overkill code to obtain different fontsize for ticks and labels
xlFontSize = get(xl,'FontSize');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 25)
set(xl, 'FontSize', xlFontSize);
yl = get(gca,'YLabel');
ylFontSize = get(yl,'FontSize');
yAY = get(gca,'YAxis');
set(yAY,'FontSize', 25)
set(yl, 'FontSize', ylFontSize);
set(gcf, 'PaperUnits', 'inches');
x_width=12.78 ;y_width=8.94;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,'Plot','epsc')

% save(savingfile,'result','S_feas','R_feas','B_feas','prange','qrange','pmin','pmax','qmin','qmax','step');

function[out] = f(x,delta,v0)
if x >= 0
    out = x/((1-delta)*v0);
else
    out = x/((1+delta)*v0);
end
end