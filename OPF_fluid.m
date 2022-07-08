function[pcon_socp,rho1,rho2] = OPF_Fluid(n,Y,Z,A,v0,delta,lambda,mu,nu)
%% SOCP
I = n;
% tic
cvx_begin sdp quiet
variables pcon_socp(1,I);
variable  W(I,I) semidefinite;
W==W'; % symmetric
W(1,1)==v0^2;

% objective function
F = cvx(zeros(1,I));
for i=2:I
%        F(i)= lambda*mu*log(pcon_socp(i));
       F(i)= lambda*mu/nu*log(pcon_socp(i)) - mu/nu*pcon_socp(i);
%        F(i)= .5*mu/lambda*pcon_socp(i)*pcon_socp(i);
end
T=sum(F);
maximize(T);

W(i,i)>=v0^2*(1-delta)^2; % volt drop
W(i,i)<=v0^2*(1+delta)^2;

WY = sum(W.*conj(Y),2); % SOCP PFE
for i=2:I
    - pcon_socp(i) == WY(i);
%     pcon_socp(i) <= lambda/mu;
end
cvx_end

z = lambda*mu*ones(1,I)./pcon_socp
% z = mu/lambda*pcon_socp.*pcon_socp
rho1= z(2);
rho2 = z(3);

end