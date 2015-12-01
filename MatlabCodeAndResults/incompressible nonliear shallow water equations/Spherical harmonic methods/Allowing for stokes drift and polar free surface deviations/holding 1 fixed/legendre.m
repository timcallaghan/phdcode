function y = legendre(x,n)
%This program is a MODIFICATION of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).

%       ===============================================
%       Purpose: Compute Legendre polynomials Pn(x)
%                and their derivatives Pn'(x)
%       Input :  x --- Argument of Pn(x)
%                n --- Degree of Pn(x) ( n = 0,1,...)
%       Output:  PN(n) --- Pn(x)
%       ===============================================

% Define the value for n=0
p0=1.0d0;
% Define the value for n=1
p1=x;

% declare pf
pf =[];

% If n is 0 or 1 we can return immediately.
if n==0
    y = p0;
elseif n==1
    y = p1;
else
    for  k=2:n;
        pf=(2.0d0.*k-1.0d0)./k.*x.*p1-(k-1.0d0)./k.*p0;
        %y=pf;
        p0=p1;
        p1=pf;
    end
    y =pf;
end
