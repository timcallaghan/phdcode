function Seta = cacheSeta(M,kappa)

%// Computes the values of the cosine basis
%// functions at each collocation point.

Seta=zeros(M,M);

deleta = pi/(M*kappa);
eps=deleta/2;
for m=1:M
	for i=1:M
		Seta(m,i)=sin((deleta*(i-1)+eps)*m*kappa);
    end
end