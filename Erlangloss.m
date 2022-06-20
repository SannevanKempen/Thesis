function[arrivalrate] = Erlangloss(n,K,lambda,nu)

rho = lambda./nu;

B = zeros(1,n);
for i=1:n
    B(i)=1;
    for k=1:K(i)
        B(i)=((rho(i)*B(i))/k)/(1+rho(i)*B(i)/k);
    end
end
arrivalrate=lambda.*(1-B);
