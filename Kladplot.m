v = 0:0.01:50;  % plotting range from -5 to 5
[x y] = meshgrid(v);  % get 2-D mesh for x and y
delta=0.1;
v0=1;
r01=0.01;
cond1 = r01/((1-delta)*v0)*x - r01/((1+delta)*v0)*y <= delta*v0;  % check conditions for these values
cond2 = -r01/((1+delta)*v0)*x + r01/((1-delta)*v0)*y <= delta*v0;
cond1 = double(cond1);  % convert to double for plotting
cond2 = double(cond2);
cond1(cond1 == 0) = NaN;  % set the 0s to NaN so they are not plotted
cond2(cond2 == 0) = NaN;
cond = cond1.*cond2;  % multiply the two condaces to keep only the common points

figure
fontsize = 35;
linewidth = 3;
surf(x,y,cond)
xlabel('p_1^+')
ylabel('p_1^-')
set(gca,'FontSize',fontsize)
view(0,90)    % change to top view
mymap1 = [0 0 1];
shading interp % turn off black edge lines
colormap(mymap1);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 20, 20], 'PaperUnits', 'Inches', 'PaperSize', [15, 10])
saveas(gca,'polyhedron_one_node','epsc') %gcf