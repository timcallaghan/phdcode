function [y,errf,errx] = mnewt(ntrial,x,tolx,tolf)
% Given and initial guess x[1..n] for a root in n dimensions, take
% ntrial Newton-Raphson steps to improve the root. Stop if the root
% converges in either summed absolute variabe increments tolx or summed
% absolute function values tolf.

% Global variables
global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp Sr Fr Ro

% Initialise the coeff vector
y=x;
errx=0;
errf=0;
n=numunknowns;
p=zeros(n,1);
delx=zeros(n,1);
fvec=zeros(n,1);
fjac=zeros(n,n);
for k=0:ntrial
	% Get the function values at y for
    % fvec and the Jacobian matrix.
    fvec=residual(y);
    fjac=fdjac(y,fvec);
    %'Condition number='
    %cond(fjac)
	%// Check function convergence
	errf=sum(abs(fvec))
	if errf <= tolf
        'Tolf reached'
        return;
    end
	%// Right hand side of linear equations.
	p = -fvec;
	% Solve the linear equations.
    delx=fjac\p;
    % Update the vector y
    y = y + delx;
    % Check for convergence
	errx=sum(abs(delx))
	if errx <= tolx
        'Tolx reached'
        return;
    end
end