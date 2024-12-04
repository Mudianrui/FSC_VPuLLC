function vir = constraintTF(e,rho)
flag = 0;%normal
if e<=-rho
    flag = -1;
elseif e>=rho
    flag = -2;
end

switch flag
    case 0
        zeta = rho^2*e/(rho+e)/(rho-e);
        vir = zeta;
    otherwise
        error(['Unhandled flag from constraint_transformation: errorflag =',num2str(flag),',e=',num2str(e),',rho=',num2str(rho),'.']);
end
