function S = cacheS(M,kappa)

%// Computes the values of the cosine basis
%// functions at each collocation point.

S=zeros(M,M);

deleta =pi/(M*kappa);
eps=deleta/2;
factor=sqrt(1.0/(2.0*pi));
for m=1:M
	for i=1:M
		S(m,i)=factor*sin((deleta*(i-1)+eps)*m*kappa);
    end
end