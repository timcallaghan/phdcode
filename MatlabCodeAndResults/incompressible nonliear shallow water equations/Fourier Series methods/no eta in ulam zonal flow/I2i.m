function I2val = I2i(i,j,a,Vzon,x,N)
%I2 Summary of this function goes here
%  Detailed explanation goes here

I2val = 0.0;
for n=1:N
    numer = 96*(-1.0+2.0*j)*(-1.0+2.0*n)*(4.0*(-1.0+j)*j+(-3.0+2.0*n)*(1.0+2.0*n));
    denom = (-3.0+2.0*j-2.0*n)*(-1.0+2.0*j-2.0*n)*(1.0+2.0*j-2.0*n)*(3.0+2.0*j-2.0*n)*(-5.0+2.0*j+2.0*n)*(-3.0+2.0*j+2.0*n)*(-1.0+2.0*j+2.0*n)*(1.0+2.0*j+2.0*n);
    I2val = I2val + x(i*N+1+n,1)*numer/denom;
end
I2val = -2.0*a*pi*I2val/Vzon;