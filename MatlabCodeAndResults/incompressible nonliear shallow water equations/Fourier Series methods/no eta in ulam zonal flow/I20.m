function I2val = I20(j,a,Vzon,x,N)
%I2 Summary of this function goes here
%  Detailed explanation goes here

I2val = 0.0;
for n=0:N
    term = ((-1)^(j-n))*(1-4*j^2-4*n^2)/(16*j^4+(1-4*n^2)^2-8*j^2*(1+4*n^2));
    I2val = I2val + x(n+1,1)*term;
end
I2val = -4.0*a*pi*I2val/Vzon;