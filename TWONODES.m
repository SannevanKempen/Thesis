function[n,v0,edges,rlist,xlist,ind_PVs] = TWONODES()
n = 3; 
v0 = 1;
ind_PVs = [2 3];

%% EDGES, RESISTANCES, REACTANCES
edges = [1 2; 2 3];
rlist = [0.01 0.01];
% xlist = zeros(1,length(edges));
xlist = 0.5*rlist;