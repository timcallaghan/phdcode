function Sphi2 = cacheSphi2(N,phi)

%// Computes the values of the cosine basis
%// functions at each collocation point.

Sphi2=zeros(N+1,N);

for n=1:N+1
	for i=1:N
		Sphi2(n,i)=sin(2*n*phi(i,1));
    end
end