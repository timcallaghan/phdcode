function y = dplgndrsimp(m,x)
	
i=0;
pmm=0;
somx2=0;
fact=0;
%// Compute dP(m,m)/dphi
pmm=x*m;
if (m>1)
	somx2=sqrt((1.0-x)*(1.0+x));
	fact=3.0;
	for i=2:m
		pmm = -fact*somx2*pmm;
		fact = fact + 2.0;
    end
end

y = pmm;
