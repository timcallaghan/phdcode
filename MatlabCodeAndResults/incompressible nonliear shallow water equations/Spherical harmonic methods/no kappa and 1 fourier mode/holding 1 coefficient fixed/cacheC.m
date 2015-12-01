function C = cacheC(M)
% Computes the values of the cosine basis
% functions at each collocation point.

% Allocate memory for the basis function storage.
C=zeros(M+1,M);

% We need M collocation points from eta=[0,2pi]
% so we have 2*pi/M=2*pi/M
deleta = pi/M;
% Use the midpoints of the grid for collocation points.
%eps=deleta/sqrt(5.0);
eps=deleta/2.0;
% The Fourier component of the Spherical Harmonics 
% nomalisation term.
normalise=sqrt(1.0/(2.0*pi));
% Populate the storage location with the basis functions
% evaluated at the collocation points.
for m=0:M
	for i=0:M-1
		C(m+1,i+1)=normalise*cos((deleta*i+eps)*m);
    end
end
