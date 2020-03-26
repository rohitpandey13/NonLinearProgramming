function f=goldarmm(func,dfdxx,dfdyy,x)
% Goldstein Armijo criteria to find the step length
con=0.0001;
dfunc=[dfdxx(x);dfdyy(x)];
d=-dfunc;
eta=0.9;
for i=1:-0.0001:0.0001
    alfa=i; % To initialize at 1
    x1=x+alfa*d;
    if (func(x1)-func(x))<= (con*alfa*d'*d)
         if (d'* [dfdxx(x1);dfdyy(x1)])  >= (eta*d'*d)
            f=i
         end
    end
end
end