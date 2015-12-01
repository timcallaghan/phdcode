function C = cacheC(M,kappa)

%// Computes the values of the cosine basis
%// functions at each collocation point.

C=zeros(M,M);

deleta =2.0*pi/(M*kappa);
eps=deleta/sqrt(5.0);
factor=sqrt(1.0/(2.0*pi));
for m=1:M
	for i=1:M
		C(m,i)=factor*cos((deleta*(i-1)+eps)*m*kappa);
    end
end
