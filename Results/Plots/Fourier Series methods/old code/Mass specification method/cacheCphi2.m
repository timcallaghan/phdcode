function Cphi2 = cacheCphi2(N,phi)

%// Computes the values of the cosine basis
%// functions at each collocation point.

Cphi2=zeros(N+1,N);

for n=1:N+1
	for i=1:N
		Cphi2(n,i)=cos(2*n*phi(i,1));
    end
end
