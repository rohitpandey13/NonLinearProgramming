%% 3. Use Cauchy's method of steepest descent to find the minimum of 
%% f(x,y) = 5x2+7y2-3xy
%% Accept a solution if gradient ration < e or if the number of 
%% iterations exceed ITMAX. Use e = 1e10-8 and ITMAX=1000
clear all
clc
syms x y
f=5*x^2+7*y^2-3*x*y;
gradient(f);
%% Given problem as a function
ff=@(x) 5*x(1)^2+7*x(2)^2-3*x(1)*x(2);
fff=@(x,y) 5*x.^2+7*y.^2-3*x.*y;
dfdx=@(x) 10*x(1) - 3*x(2);
dfdy=@(x) 14*x(2) - 3*x(1);
%% Lambda calculation using Goldstein Armijo criteria
%% Xk+1 = Xk + lambda * dk
%% Choosing lambda > 0  
% Following criteria's are used to check the lambda > 0
% 1. f(Xk+1)-f(Xk) <= alpha*lambda * dk' * gradf(Xk)
% 2. dk'*gradf(Xk+1) >= eta * dk' * gradf(Xk)
%% Let initial lambda=0.1
%% Let initial guess x0 = [0;-1]
fcontour(fff,[-4 4 -2 10],'LevelStep',10)
grid
x0=0;y0=-1;
hold on
plot(x0,y0,'rx');
x0 =[0;-1];
d0=-[dfdx(x0);dfdy(x0)];
lam0=0.1;
x1=x0+lam0*d0;
i=0;
ITMAX=1000;
betaa1=100;
epsilon=0.00000001;
tic
while (i < ITMAX) && (betaa1 >= epsilon)
    %lam1=goldarmm(ff,dfdx,dfdy,x1);
    dfunc=[dfdx(x1);dfdy(x1)];
    lam1=linesearch(ff,dfunc,x1);
    x2=x1+lam1*-[dfdx(x1);dfdy(x1)];
    plot(x2(1,1),x2(2,1),'rx');
    pause(0.05);
    betaa1=((dfdx(x2)).^2+(dfdy(x2)).^2)/(1+abs(ff(x2)));
    x0=x1;
    x1=x2;
    i=i+1;
end
toc
hold off
disp("Local minimum solution approach: Cauchy's steepest descent method")
fprintf("The solution is obtained after %d iterations\n",i)
fprintf("The objective function value is %d\n",ff(x2))
fprintf("The local minimum point : %0.5f\n",x1)