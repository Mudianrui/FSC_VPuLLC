function [mu,upsilon] = constraintTF_pd(e,rho)
flag = 0;%normal
if e<=-rho
    flag = -1;
elseif e>=rho
    flag = -2;
end

switch flag
    case 0
%         zeta = rho^2*e/(rho+e)/(rho-e);
        mu = rho^2*(rho^2+e^2)/((rho+e)*(rho-e))^2;
        upsilon = -2*rho*e^3/((rho+e)*(rho-e))^2;
    otherwise
        error(['Unhandled flag from constraint_transformation: errorflag =',num2str(flag),',e=',num2str(e),',rho=',num2str(rho),'.']);
end
