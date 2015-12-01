function C = cacheC(M)

%// Computes the values of the cosine basis
%// functions at each collocation point.

C=zeros(M+1,M);

deleta = pi/M;
eps=deleta/2.0;
factor=sqrt(1.0/(2.0*pi));
for m=0:M
	for i=0:M-1
		C(m+1,i+1)=factor*cos((deleta*i+eps)*m);
    end
end
