function P = cacheL(mu,M,kappa)
% Computes the values of the associated Legendre
% functions at each collocation point. It also includes
% the normalising term for the spherical harmonics 
% to reduce the size of the function values....see plgndr for
% the normalising term.

% Allocate memory for the basis function storage.
P=zeros(M+2,2*M,M);

% Populate the storage location with the basis functions
% evaluated at the collocation points. First assign all 
% those basis functions that are NOT using basis recombination
for m=0:M
	for l=m:2*M+m-1
		for k=1:M
			P(m+1,l-m+1,k)=plgndr(l,m,mu(k,1));
        end
    end
end

% Now assign those basis functions that ARE using
% basis recombination
normalize=0.0;
for l=2:2*M+1
	for k=1:M
        normalize = sqrt((2.0*l+1.0)/2.0);
        % Check to see if l is odd or even
        if (mod(l,2)==0)
            % Since l is even we subtract off the
            % even basis function (which is just unity)
		    P(M+2,l-1,k)=normalize*(legendrep(mu(k,1),l)-1.0);
        else
            % l must be odd so we now subtract off the
            % first odd basis function (which is mu=sin(phi)
            P(M+2,l-1,k)=normalize*(legendrep(mu(k,1),l)-mu(k,1));
        end
    end
end
