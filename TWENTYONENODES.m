function[n,v0,edges,r,x,ind_PVs] = TWENTYONENODES()
n = 22;
v0 = 1;
ind_PVs = 2:22;
edges = [ 1     2
     2     3
     2     4
     4     5
     4     9
     5     6
     6     7
     6     8
     9    10
     9    11
    11    12
    11    13
    13    14
    14    15
    14    16
    16    17
    17    18
    17    19
    19    20
    20    21
    20    22];


r = [0.0030    0.0005    0.0045    0.0016    0.0061    0.0108    0.0005    0.0024    0.0005    0.0056    0.0005    0.0033    0.0086    0.0002    0.0005    0.0027    0.0008    0.0047    0.0011    0.0007    0.0044];
x = [0.0015    0.0002    0.0023    0.0008    0.0032    0.0056    0.0003    0.0012    0.0002    0.0029    0.0002    0.0017    0.0045    0.0001    0.0002    0.0014    0.0004    0.0024    0.0005    0.0004    0.0023];
end