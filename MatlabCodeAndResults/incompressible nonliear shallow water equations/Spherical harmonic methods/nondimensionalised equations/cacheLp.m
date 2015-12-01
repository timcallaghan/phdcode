function Pp=cacheLp(phi,M,kappa)

Pp=zeros(M,2*M,M);


%// Computes the derivatives of the associated Legendre
%// functions at each collocation point. It also includes
%// the normalising term for the sperical harmonics 
%// to reduce the size of the function values.
m=0;
l=0;
k=0;

for m=1:M
	for l=kappa*m:2*M+kappa*m-1
		for k=1:M
			Pp(m,l-m*kappa+1,k)=dplgndr(l,m*kappa,phi(k,1));
        end
    end
end
