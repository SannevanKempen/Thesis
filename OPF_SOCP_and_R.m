function[pcon_mod,qcon_mod,pgen_mod,qgen_mod,pcon_socp,qcon_socp,pgen_socp,qgen_socp,volt_socp,W,WY,time_mod,time_socp] = OPF_SOCP_and_R(n,Y,A,v0,delta,pconub,qconub,pgenub,qgenub);
%% modified
% tic
cvx_begin sdp quiet
variables pcon_mod(1,n);
variables qcon_mod(1,n);
variables pgen_mod(1,n);
variables qgen_mod(1,n);

% objective function
F = cvx(zeros(1,n));
for i=1:n
        F(i)= log(pcon_mod(i)); 
end
T=sum(F);
maximize(T);


%constraints
A*[pcon_mod qcon_mod pgen_mod qgen_mod]' <= delta*v0
eye(4*n)*[pcon_mod qcon_mod pgen_mod qgen_mod]' >= 0;
eye(n)*pcon_mod' <= pconub';
eye(n)*qcon_mod' <= qconub';
eye(n)*pgen_mod' <= pgenub';
eye(n)*qgen_mod' <= qgenub';
cvx_end
time_mod = toc;

%% SOCP
I = n+1;
pconub = [10^3 pconub];
qconub = [10^3 qconub];
pgenub = [10^3 pgenub];
qgenub = [10^3 qgenub];

% tic
cvx_begin sdp quiet
variables pcon_socp(1,I);
variables pgen_socp(1,I);
variables qcon_socp(1,I);
variables qgen_socp(1,I);
variable  W(I,I) semidefinite;
W==W'; % symmetric
W(1,1)==v0^2;

% objective function
F = cvx(zeros(1,I));
for i=2:I
        F(i)= log(pcon_socp(i)); 
end
T=sum(F);
maximize(T);

% constraints
eye(4*I)*[pcon_socp qcon_socp pgen_socp qgen_socp]' >= 0; % admissible set
eye(I)*pcon_socp' <= pconub';
eye(I)*qcon_socp' <= qconub';
eye(I)*pgen_socp' <= pgenub';
eye(I)*qgen_socp' <= qgenub';

W(i,i)>=v0^2*(1-delta)^2; % volt drop
W(i,i)<=v0^2*(1+delta)^2;

WY = sum(W.*conj(Y),2); % SOCP PFE
for i=2:I
    pgen_socp(i) - pcon_socp(i) + 1i*(qgen_socp(i)-qcon_socp(i)) == WY(i);
end
cvx_end
time_socp = toc;

volt_socp = diag(W);


end