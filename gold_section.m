function f=gold_section(func,low,high,tol)
tol=Inf;
flow=func(low);
fhigh=func(high);
r=0.5*(sqrt(5)-1);
T=1-r;
Lorg=high-low;
x1=(1-T)*low+T*high;
x2=T*low+(1-T)*high;
iter=0;
tic
while tol > 0.00000001
    if func(x2) > func(x1)
        high=x2;
        x2=x1;
        x1=(1-T)*low+T*high;
        tol= (high-low) / Lorg;
    else
        low=x1;
        x1=x2;
        x2=T*low+(1-T)*high;
        tol = (high-low)/Lorg;
    end
    iter=iter+1;
end
toc
result=(low+high)/2;
disp("Local minimum solution approach: Golden Section method")
fprintf("The solution is obtained after %d iterations\n",iter)
fprintf("The local minimum point using Algorithm 2 : %0.5f\n",result) 