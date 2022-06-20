function [pcon,pgen,W] = OPF_SOCP(n,Y,v0,delta,pconub,pgenub)

I = n+1;
pconub = [10^3 pconub];
pgenub = [10^3 pgenub];

cvx_begin sdp quiet
variables pcon(1,I);
variables pgen(1,I);
variable  W(I,I) semidefinite;
W==W'; % symmetric
W(1,1)==v0^2;

%objective function
F=cvx(zeros(1,I));
for i=2:I
    if pconub(i) > 0 % dont want to take log of 0
        F(i)=log(pcon(i));
        %         F(i) = pcon(i);
    end
end
T=sum(F);
maximize(T);

% constraints
eye(2*I)*[pcon pgen]' >= 0; % admissible set
eye(I)*pcon' <= pconub';
eye(I)*pgen' <= pgenub';

W(i,i)>=v0^2*(1-delta)^2; % volt drop 
W(i,i)<=v0^2*(1+delta)^2;

YW = sum(Y.*W,2); % SOCP PFE
for i=1:I
    pgen(i) - pcon(i) == YW(i);
end

cvx_end
end