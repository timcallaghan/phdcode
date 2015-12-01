function Cphi1 = cacheCphi1(N,phi)

%// Computes the values of the cosine basis
%// functions at each collocation point.

Cphi1=zeros(N+1,N);

for n=1:N+1
	for i=1:N
		Cphi1(n,i)=cos((2*n-1)*phi(i,1));
    end
end
