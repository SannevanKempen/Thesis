function[pcon_mod,pgen_mod,qgen_mod,pcon_socp,pgen_socp,qgen_socp,volt_socp,YW,W,g] = Test(n,Y,A,v0,delta,pconub,pgenub,qgenub,lambda,mu,nu,K,c);
%% modified

lambdaK = Erlangloss(n,K,lambda,nu);
f = @(x) lambdaK.*x./(x.*mu+nu); % Markovian model!
g = f(c);

cvx_begin sdp quiet
variables pcon_mod(1,n);

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
A*[pcon_mod zeros(1,n) zeros(1,n) zeros(1,n)]' <= delta*v0
eye(n)*[pcon_mod]' >= 0;
eye(n)*pcon_mod' <= pconub';
eye(n)*pcon_mod' <= g';
cvx_end

%% SOCP
I = n+1;
pconub = [10^3 pconub];
lambda = [0 lambda];
mu = [0 mu];
nu = [0 nu];
K = [0 K];

lambdaK = Erlangloss(I,K,lambda,nu);
f = @(x) lambdaK.*x./(x.*mu+nu); % Markovian model!
g = f(c);

cvx_begin sdp quiet
variables pcon_socp(1,I);
variable  W(I,I);
W >= 0; % semi-definite
W==W'; % symmetric
W(1,1)==v0^2;

% objective function
F = cvx(zeros(1,I));
% w(1) = 0;
% w(2) = 0.01;
% w(3) = 0.015;
for i=2:I
    if lambda(i)>0
%         F(i)= (lambdaK(i)/nu(i))*w(i)*log(pcon_socp(i))-mu(i)*w(i)*pcon_socp(i)/nu(i);
        F(i)= (lambdaK(i)/nu(i))*log(pcon_socp(i))-mu(i)*pcon_socp(i)/nu(i); % Markovian model!
    end
end
T=sum(F);
maximize(T);

% constraints
eye(I)*[pcon_socp]' >= 0; % admissible set
eye(I)*pcon_socp' <= pconub';

W(i,i)>=v0^2*(1-delta)^2; % volt drop
W(i,i)<=v0^2*(1+delta)^2;

YW = sum(W.*conj(Y),2); % SOCP PFE
for i=2:I
    -pcon_socp(i) == YW(i);
    if lambda(i) > 0
        pcon_socp(i) <= g(i);
    else
        pcon_socp(i) == 0
    end
end
cvx_end

volt_socp = diag(W);
pgen_mod = zeros(1,n);
pgen_socp = zeros(1,I);
qgen_mod = zeros(1,n);
qgen_socp = zeros(1,I);
end