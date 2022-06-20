function [pcon,pgen] = OPF_modified(n,A,v0,delta,pconub,pgenub)
% numbering: substation is node 0

%Optimazation problem
cvx_begin sdp quiet
% control variables
variables pcon(1,n);
variables pgen(1,n);

%objective function
F=cvx(zeros(1,n));
%
for i=1:n
    if pconub(i) > 0
        F(i)=log(pcon(i));
%         F(i) = pcon(i);
    end
end
T=sum(F);
maximize(T);

%constraints
A*[pcon pgen]' <= delta*v0
eye(2*n)*[pcon pgen]' >= 0;
eye(n)*pcon' <= pconub';
eye(n)*pgen' <= pgenub';
cvx_end

end