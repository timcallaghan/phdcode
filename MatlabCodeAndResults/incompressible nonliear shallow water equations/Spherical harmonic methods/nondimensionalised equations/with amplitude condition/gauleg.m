function mu = gauleg(x1,x2,M)

%/* Given the lower and upper limits of integration x1 and
% x2, this routine returns the array x[0..n-1] of
% length n, containing the abscissas of the 
% Gauss-Legendre n-point quadrature formula.
% */
EPS=1.0e-14;  %// relative precision
% Variables used in the calculations
m=0;
j=0;
i=0;
z1=0;
z=0;
xm=0;
xl=0;
pp=0;
p3=0;
p2=0;
p1=0;
mu=zeros(M,1);

n=2*M;
%// The roots are symmetric in the interval
%// so we need only find half of them. However
%// we require M roots in half the interval so
%// we just compute these...
m=(n+1)/2;
xm=0.5*(x2+x1);
xl=0.5*(x2-x1);
for i=0:m-1
	%// Loop over the desired roots.
	z=cos(pi*(i+0.75)/(n+0.5));
	%// Starting with this approximation to
	%// the root, we proceed with Newton's
	%// method and refine this guess.
	
    % Should have a do/while loop here but it seems MATLAB
    % doesn't account for this control structure...???
    % So just do first part of do loop here
    p1=1.0;
	p2=0.0;
	for j=0:n-1
		%// Loop up the recurrence relation to get
		%// the Legendre polynomial evaluated at z.
		p3=p2;
		p2=p1;
		p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
    end
	%// p1 is now the desired Legendre polynomial, Now compute
	%// pp, it's derivative by a standard relation involving
	%// also p2, the polynomial of lower order.
	pp=n*(z*p1-p2)/(z*z-1.0);
	z1=z;
	%// Newton's method...
	z=z1-p1/pp;
    % Now perform the while part of the loop
    while (abs(z-z1) > EPS)
		p1=1.0;
		p2=0.0;
		for j=0:n-1
			%// Loop up the recurrence relation to get
			%// the Legendre polynomial evaluated at z.
			p3=p2;
			p2=p1;
			p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
        end
		%// p1 is now the desired Legendre polynomial, Now compute
		%// pp, it's derivative by a standard relation involving
		%// also p2, the polynomial of lower order.
		pp=n*(z*p1-p2)/(z*z-1.0);
		z1=z;
		%// Newton's method...
		z=z1-p1/pp;
    end    
	%// Scale the root to the desired interval.
	%//x[i+1]=xm-xl*z;
	%// Put in it's symmetric counterpart.
	mu(n/2-i,1)=xm+xl*z;
end
