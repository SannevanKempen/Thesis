eta = .9;
gamma = .1;
delta = .1;

n = 2; % 3;
v0 = 1;
edges = zeros(n-1,2);
edges(1,1) = 1; % [1 2; 2 3];
edges(1,2) = 2;
r = [0.05]; % [.05 .05];
x = eta*r;
ind_PVs = [];
% matpower_name = 'matpower_cases/case_two_nodes.m';
matpower_name = 'matpower_cases/case_one_node.m';
p = 1;
q = gamma*p;
s = p+1i*q;
Ymat = zeros(n,n);
for i=1:n-1
    Ymat(edges(i,1),edges(i,2)) = -1/(r(i)+1i*x(i));
    Ymat(edges(i,2),edges(i,1)) = -1/(r(i)+1i*x(i));
end
for i=1:n
    Ymat(i,i) = -sum(Ymat(i,:));
end
Zmat = inv(Ymat(2:end,2:end));
Rmat = real(Zmat);
Xmat = imag(Zmat);
Gmat = real(Ymat);
Bmat = imag(Ymat);
mpc = loadcase(matpower_name);
mpc.gen(1, 6) = v0; % set correct values for v0,r and x in mpc
for j=1:n-1
    mpc.branch(j, 3) = r(j);
    mpc.branch(j, 4) = x(j);
    mpc.bus(j+1,3) = p;
    mpc.bus(j+1,4) = q;
end

mag = zeros(1,n);
angles = zeros(1,n);
[results,succes] = runpf(mpc,mpoption('verbose',0,'out.all',0));
if succes == 1
    nr_of_violations = sum(results.bus(:,8) < (1-delta)*v0) + sum(results.bus(:,8) > (1+delta)*v0);
    if nr_of_violations == 0
        mag = results.bus(:,8)'
        angles = deg2rad(results.bus(:,9))'
        
        cosangles = cos(angles)
        sinangles = sin(angles)
        
        pc = max(p,0);
        pg = max(-p,0);
        
        cosbounds = [((1-delta)^2*v0^2+(1-eta*gamma)*r*pc-r*pg)/((1+delta)*v0^2) ((1+delta)^2*v0^2+r*pc+(-1+eta*gamma)*r*pg)/((1-delta)*v0^2)]
        sinbounds = [(-gamma*r*pg - eta*r*pc)/((1+delta)*v0^2) (gamma*r*pc+eta*r*pg)/((1-delta)*v0^2)]
    end
end




% sind(meshgrid(angles) - meshgrid(angles)');
% cosd(meshgrid(angles) - meshgrid(angles)');
%
% V0 = mag(1)*exp(1i*(angles(1)));
% V1 = mag(2)*exp(1i*(angles(2)));
% V2 = mag(3)*exp(1i*(angles(3)));
%
% V1*conj(V0)*conj(Ymat(2,1))+V1*conj(V1)*conj(Ymat(2,2))+V1*conj(V2)*conj(Ymat(2,3));
% V2*conj(V0)*conj(Ymat(3,1))+V2*conj(V1)*conj(Ymat(3,2))+V2*conj(V2)*conj(Ymat(3,3));
%
%
% mag(1)*mag(2)*(Bmat(2,1)*sin(angles(2)-angles(1))+Gmat(2,1)*cos(angles(2)-angles(1)))+mag(2)*mag(2)*(Bmat(2,2)*sin(angles(2)-angles(2))+Gmat(2,2)*cos(angles(2)-angles(2))) +mag(3)*mag(2)*(Bmat(2,3)*sin(angles(2)-angles(3))+Gmat(2,3)*cos(angles(2)-angles(3)));
% sin(angles(2)-angles(1));;
%
% (mag(1)*Gmat(2,1)*(sin(angles(1))+eta*cos(angles(1))) + mag(2)*Gmat(2,2)*(sin(angles(2))+eta*cos(angles(2))) + mag(3)*Gmat(2,3)*(sin(angles(3))+eta*cos(angles(3))))/(mag(1)*Gmat(2,1)*(-eta*sin(angles(1))+cos(angles(1))) + mag(2)*Gmat(2,2)*(-eta*sin(angles(2))+cos(angles(2))) + mag(3)*Gmat(2,3)*(-eta*sin(angles(3))+cos(angles(3))));
% (mag(1)*Gmat(3,1)*(sin(angles(1))+eta*cos(angles(1))) + mag(2)*Gmat(3,2)*(sin(angles(2))+eta*cos(angles(2))) + mag(3)*Gmat(3,3)*(sin(angles(3))+eta*cos(angles(3))))/(mag(1)*Gmat(2,1)*(-eta*sin(angles(1))+cos(angles(1))) + mag(2)*Gmat(2,2)*(-eta*sin(angles(2))+cos(angles(2))) + mag(3)*Gmat(2,3)*(-eta*sin(angles(3))+cos(angles(3))));
%












