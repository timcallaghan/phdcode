function Sphi1 = cacheSphi1(N,phi)

%// Computes the values of the cosine basis
%// functions at each collocation point.

Sphi1=zeros(N+1,N);

for n=1:N+1
	for i=1:N
		Sphi1(n,i)=sin((2*n-1)*phi(i,1));
    end
end