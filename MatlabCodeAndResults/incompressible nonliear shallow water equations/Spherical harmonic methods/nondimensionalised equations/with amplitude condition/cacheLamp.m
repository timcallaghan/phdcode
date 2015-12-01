function Pamp = cacheLamp(M,kappa)

%// Computes the values of the associated Legendre
%// functions at each collocation point. It also includes
%// the normalising term for the sperical harmonics 
%// to reduce the size of the function values....see plgndr for
%// the normalising term.
int m0;
l=0;

Pamp=zeros(M,2*M);

%// The xvalue at lat pi/4 ie x=sin(pi/4)
xval=sin(pi/4.0);
	
for m=1:M
	for l=kappa*m:2*M+kappa*m-1
        Pamp(m,l-m*kappa+1)=plgndr(l,m*kappa,xval);
    end
end