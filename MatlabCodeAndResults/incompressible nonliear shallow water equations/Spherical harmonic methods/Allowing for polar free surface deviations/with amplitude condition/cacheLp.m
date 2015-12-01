function Pp=cacheLp(phi,M,kappa)

%// Computes the derivatives of the associated Legendre
%// functions at each collocation point. It also includes
%// the normalising term for the spherical harmonics 
%// to reduce the size of the function values.

Pp=zeros(M+1,2*M,M);


for m=0:M
	for l=kappa*m:2*M+kappa*m-1
		for k=1:M
			Pp(m+1,l-m*kappa+1,k)=dplgndr(l,m*kappa,phi(k,1));
        end
    end
end
