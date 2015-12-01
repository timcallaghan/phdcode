function P = cacheL(mu,M,kappa)

%// Computes the values of the associated Legendre
%// functions at each collocation point. It also includes
%// the normalising term for the spherical harmonics 
%// to reduce the size of the function values....see plgndr for
%// the normalising term.

P=zeros(M,2*M,M);

for m=1:M
	for l=kappa*m:2*M+kappa*m-1
		for k=1:M
			P(m,l-m*kappa+1,k)=plgndr(l,m*kappa,mu(k,1));
        end
    end
end