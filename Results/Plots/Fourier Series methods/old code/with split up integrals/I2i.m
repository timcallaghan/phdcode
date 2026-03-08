function I2val = I2i(i,j,a,Vzon,x,N)
%I2 Summary of this function goes here
%  Detailed explanation goes here

I2val = 0.0;
for n=1:N
    numer = 96*(-1+2*j)*(-1+2*n)*(4*(-1+j)*j+(-3+2*n)*(1+2*n));
    denom = (-3+2*j-2*n)*(-1+2*j-2*n)*(1+2*j-2*n)*(3+2*j-2*j)*(-5+2*j+2*n)*(-3+2*j+2*n)*(-1+2*j+2*n)*(1+2*j+2*n);
    I2val = I2val + x(i*N+1+n,1)*numer/denom;
end
I2val = -2.0*a*pi*I2val/Vzon;