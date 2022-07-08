function[alc_mod,alc_socp,alc_ldf] = OPF_SOCP_and_R_and_LDF(n,Y,Z,A,v0,delta,u,P2gen,Q2gen,z)
%% modified
% tic
cvx_begin sdp quiet
variables alc_mod(1,n);

% objective function
F = cvx(zeros(1,n));
for i=1:n
        F(i)= u(i)*z(i)*log(alc_mod(i));
end
T=sum(F);
maximize(T);
pcon_mod = z.*alc_mod;
qcon_mod = zeros(1,n);
pgen_mod = zeros(1,n);
qgen_mod = zeros(1,n);
pgen_mod(end) = P2gen;
qgen_mod(end) = Q2gen;

%constraints
A*[pcon_mod qcon_mod pgen_mod qgen_mod]' <= delta*v0
cvx_end


%% LinDistFlow
cvx_begin sdp quiet
variables alc_ldf(1,n);
variable  Vlinsq(1,n);

% objective function
F = cvx(zeros(1,n));
for i=1:n
         F(i)= u(i)*z(i)*log(alc_ldf(i));
end
T=sum(F);
maximize(T);

pcon_ldf = z.*alc_ldf;
qcon_ldf = zeros(1,n);
pgen_ldf = zeros(1,n);
qgen_ldf = zeros(1,n);
pgen_ldf(end) = P2gen;
qgen_ldf(end) = Q2gen;

Vlinsq(i)>=v0^2*(1-delta)^2; % volt drop
Vlinsq(i)<=v0^2*(1+delta)^2;

Zs = real(Z)*[pcon_ldf - pgen_ldf]' + imag(Z)*[qcon_ldf - qgen_ldf]';
for i=1:n
    Vlinsq(i) == v0^2 - 2*Zs(i);
end
cvx_end


%% SOCP
I = n+1;
u = [0 u];
z = [0 z];

% tic
cvx_begin sdp quiet
variables alc_socp(1,I);
variable  W(I,I) semidefinite;
W==W'; % symmetric
W(1,1)==v0^2;

pcon_socp = z.*alc_socp;
qcon_socp = zeros(1,I);
pgen_socp = zeros(1,I);
qgen_socp = zeros(1,I);
pgen_socp(end) = P2gen;
qgen_socp(end) = Q2gen;


% objective function
F = cvx(zeros(1,I));
for i=2:I
        F(i)= u(i)*z(i)*log(alc_socp(i));
end
T=sum(F);
maximize(T);

W(i,i)>=v0^2*(1-delta)^2; % volt drop
W(i,i)<=v0^2*(1+delta)^2;

WY = sum(W.*conj(Y),2); % SOCP PFE
for i=2:I
    pgen_socp(i) - pcon_socp(i) + 1i*(qgen_socp(i)-qcon_socp(i)) == WY(i);
end
cvx_end


end