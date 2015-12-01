function P = cacheL(mu,M,kappa)

%// Computes the values of the associated Legendre
%// functions at each collocation point. It also includes
%// the normalising term for the spherical harmonics 
%// to reduce the size of the function values....see plgndr for
%// the normalising term.

P=zeros(M+2,2*M,M);

% Assign all those basis functions that are not
% using recombination
for m=0:M
	for l=kappa*m:2*M+kappa*m-1
		for k=1:M
			P(m+1,l-m*kappa+1,k)=plgndr(l,m*kappa,mu(k,1));
        end
    end
end

% Now assign those basis functions that are using
% basis recombination
normalize=0.0;
for l=2:2*M+1
	for k=1:M
        normalize = sqrt((2.0*l+1.0)/2.0);
		P(M+2,l-1,k)=normalize*(legendre(mu(k,1),l)-1.0);
    end
end
