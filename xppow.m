function ff=xppow(x)
rb = 0.4;
r1 = 2;
rt = 3;
b = (r1-rb)/(rt-r1);
a = (rt-rb)*b;
p = a*x^2/(1+b*x^2)+rb;
ff = sgn(x,p);
