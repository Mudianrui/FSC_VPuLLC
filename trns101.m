function [s] = trns101(x,k)
if abs(x)<=1
    s = 0;
else
    f = k/(1-x^2);
    s = exp(f);
end
