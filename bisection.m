function f = bisection(funct,low,high,tol)
tic,syms f(x);
f(x)=funct;
pretty(f)
deriv=diff(f);
err=inf;
itr=0;
flow=deriv(low);
fhigh=deriv(high);
if flow < 0 && fhigh > 0
    while err > tol
       m=(low+high)/2;
       fdmid=deriv(m);
        if fdmid < 0
             low=m;
        elseif fdmid>0
            high=m;
        else
            xopt=m;
        end
    itr=itr+1;
    err=abs(high-low);
    end
else
    disp("Find a new low and high values")
    return
end
xopt=m;toc
disp("Local minimum solution approach: Golden Bisection method")
fprintf("The solution is obtained after %d iterations\n",itr)
fprintf("The local minimum point : %0.5f\n",xopt)
    
