function T1=smoothTfun1(z)
if z>0
    T1=exp((z-1)/z);
else
    T1=0;
end
