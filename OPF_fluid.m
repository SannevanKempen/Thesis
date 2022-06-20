function[pcon_mod,pgen_mod,qgen_mod,pcon_socp,pgen_socp,qgen_socp,volt_socp] = OPF_fluid(n,Y,A,v0,delta,pconub,pgenub,qgenub,lambda,mu,nu,K,c);
%% modified
lambdaK = Erlangloss(n,K,lambda,nu);
f = @(x) lambdaK.*x./(x.*mu+nu); % Markovian model!
g = f(c);

cvx_begin sdp quiet
variables pcon_mod(1,n);
variables pgen_mod(1,n);
variables qgen_mod(1,n);

% objective function
F = cvx(zeros(1,n));
for i=1:n
    if lambda(i)>0
        F(i)= (lambdaK(i)/nu(i))*log(pcon_mod(i))-mu(i)*pcon_mod(i)/nu(i); % Markovian model!
    end
end
T=sum(F);
maximize(T);

%constraints
A*[pcon_mod zeros(1,n) pgen_mod qgen_mod]' <= delta*v0
eye(3*n)*[pcon_mod pgen_mod qgen_mod]' >= 0;
eye(n)*pcon_mod' <= pconub';
eye(n)*pgen_mod' <= pgenub';
eye(n)*qgen_mod' <= qgenub';
eye(n)*pcon_mod' <= g';
cvx_end

%% SOCP
I = n+1;
pconub = [10^3 pconub];
pgenub = [10^3 pgenub];
qgenub = [10^3 qgenub];
lambda = [0 lambda];
mu = [0 mu];
nu = [0 nu];
K = [0 K];

lambdaK = Erlangloss(I,K,lambda,nu);
f = @(x) lambdaK.*x./(x.*mu+nu); % Markovian model!
g = f(c);

cvx_begin sdp quiet
variables pcon_socp(1,I);
variables pgen_socp(1,I);
variables qgen_socp(1,I);
variable  W(I,I) semidefinite;
W==W'; % symmetric
W(1,1)==v0^2;

% objective function
F = cvx(zeros(1,I));
for i=2:I
    if lambda(i)>0
        F(i)= (lambdaK(i)/nu(i))*log(pcon_socp(i))-mu(i)*pcon_socp(i)/nu(i); % Markovian model!
    end
end
T=sum(F);
maximize(T);

% constraints
eye(3*I)*[pcon_socp pgen_socp qgen_socp]' >= 0; % admissible set
eye(I)*pcon_socp' <= pconub';
eye(I)*pgen_socp' <= pgenub';
eye(I)*qgen_socp' <= qgenub';

W(i,i)>=v0^2*(1-delta)^2; % volt drop
W(i,i)<=v0^2*(1+delta)^2;

YW = sum(conj(Y).*W,2); % SOCP PFE
for i=1:I
    pgen_socp(i) - pcon_socp(i) + 1i*qgen_socp(i) == YW(i)
    if lambda(i) > 0
        pcon_socp(i) <= g(i);
    end
end
cvx_end

volt_socp = diag(W);

end