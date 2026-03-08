function y = basepmass(gamma,a,Ma,Fr)
% Computes the base mass for the compressible
% shallow water equations with a constant free
% surface height.

B=(gamma-1)/gamma*(Ma/Fr)^2;
f1=B*(a+1);

y=2*(f1*(gamma-1)/a)^2+2*B*f1*(gamma-1)*gamma/a+gamma*(2*gamma-1)*B^2;

