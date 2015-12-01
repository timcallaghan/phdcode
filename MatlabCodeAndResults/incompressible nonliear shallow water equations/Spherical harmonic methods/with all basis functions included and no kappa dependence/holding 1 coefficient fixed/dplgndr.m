function y = dplgndr(l,m,phival)
% ********NOTE*********
% This function accepts phi NOT x as argument!
% *********************
% Computes the associated Legrendre polynomials
% derivatives dP(m,l)/dphi. Here m and l are integers satisfying
% 0<=m<=l, while phi lies in the range -pi/2<=x<=pi/2

normalize = sqrt((2.0*l+1.0)/2.0*factrl(l-m)/factrl(l+m));

% Check if m==0, if so we just have the Legendre polynomials
if (m==0)
    normalize = sqrt((2.0*l+1.0)/2.0);
    y=normalize*dlegendre(l,phival);
    return;
else
    %// If l==m then we just use the simple formula for
    %// the derivative...
    if (l==m)
	    y = normalize*dplgndrsimp(m,sin(phival));
        return;
    else
	    %// else we don't have l==m so we use the more
	    %// general formula for the calculation
	    y = normalize*((l+m)*plgndrmod(l-1,m,sin(phival))-l*sin(phival)*plgndrmod(l,m,sin(phival)));
        return;
    end
end