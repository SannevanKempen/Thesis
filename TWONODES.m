function[n,v0,edges,rlist,xlist,ind_PVs] = TWONODES()
n = 3; 
v0 = 1;
ind_PVs = [3];

%% EDGES, RESISTANCES, REACTANCES
edges = [1 2; 2 3];
rlist = [0.01 0.005];
xlist = zeros(1,length(edges));
% xlist = rlist;