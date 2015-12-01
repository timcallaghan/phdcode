function Pamp = cacheLamp(M,kappa)

%// Computes the values of the associated Legendre
%// functions at each collocation point. It also includes
%// the normalising term for the sperical harmonics 
%// to reduce the size of the function values....see plgndr for
%// the normalising term.

Pamp=zeros(M+1,2*M);

%// The xvalue at lat pi/2 ie x=sin(pi/2)=1
xval=1.0;
	
for m=0:M
	for l=kappa*m:2*M+kappa*m-1
        Pamp(m+1,l-m*kappa+1)=plgndr(l,m*kappa,xval);
    end
end