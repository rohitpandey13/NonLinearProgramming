%% OR 506
%% Homework 3
%% Submitted by Srinivasan Balan (200209659)
%% 1. Find the diagonal matrix E so that A+E=LDLT, D a diagonal matrix 
%% with positive diagonal entries, L a lower triangular matrix and 
clear all
clc
A = [1 4 3;4 2 5;3 5 3]
%% Refer handout for the solution
E=[1.41 0 0;0 4.66 0;0 0 0.76]
D=[1.25 0 0; 0 0.38 0; 0 0 0.267]
% Used the Levenberg Marquadt solution approach to calculate r1,r2,r3
%% 2. Let f(x,y)=2x2+y2-2xy+2x3+x4
%% Minimize the function starting from x0 = (0,-2)T.
%% Verify that d0 = (0,1)T is a direction of descent
clear all
clc
syms x y
f=2*x^2+y^2-2*x*y+2*x^3+x^4
gradient(f)
%% Given problem as a function
ff=@(x,y) 2*x.^2+y.^2-2*x.*y+2*x.^3+x.^4;
dfdx=@(x,y) 4*x^3 + 6*x^2 + 4*x - 2*y;
dfdy=@(x,y) 2*y - 2*x;
%% Direction of descent dk using Goldstein Armijo criteria
%% Xk+1 = Xk + lambda * dk
%% Choosing lambda > 0  
% Following criteria's are used to check the lambda > 0
% 1. f(Xk+1)-f(Xk) <= alpha*lambda * dk' * gradf(Xk)
% 2. dk'*gradf(Xk+1) >= eta * dk' * gradf(Xk)
%% Let initial lambda=0.5
fcontour(ff,[-4 4 -2 10],'LevelStep',10)
grid
x0=0;y0=-2;
hold on
plot(x0,y0,'ro');
x0 =[0;-2];
d0=-[dfdx(0,-2);dfdy(0,-2)];
lam0=0.1;
x1=x0+lam0*d0;
alpha=0.0001;
eta=0.9;
i=0;
ITMAX=100;
beta1=100;
beta2=100;
epsilon=0.001;
lam1=0.1;
tic
while i < ITMAX && (beta1 > epsilon || beta2 > epsilon)
    lam1=0.1+(0.1)*rand(1,1);
    if (ff(x1(1,1),x1(2,1))-ff(x0(1,1),x0(2,1))) <= (alpha*lam1*-[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))]'*[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))])
        if -[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))]'*[dfdx(x1(1,1),x1(2,1));dfdy(x1(1,1),x1(2,1))] >= (eta*-[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))]'*[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))])
            x2=x1+lam1*-[dfdx(x1(1,1),x1(2,1));dfdy(x1(1,1),x1(2,1))];
            plot(x2(1,1),x2(2,1),'ro');
            pause(0.05);
            beta1=abs(x2(1,1)-x1(1,1));
            beta2=abs(x2(2,1)-x2(2,1));
            x0=x1;
            x1=x2;
            i=i+1;
        end
    end
end
toc
hold off
disp("Step length calculation using: Goldstein Armijo method")
disp("Local minimum solution approach: Steepest descent method")
fprintf("The solution is obtained after %d iterations\n",i)
fprintf("The local minimum point : %0.5f\n",x1)
%% Verifying d0 = (0,1)T as a direction descent
a=dfdx(0,1)
b=dfdy(0,1)
gradf=[a;b]
dk=-gradf
dd=dk'*gradf
if dd < 0
    disp("The dk = (0,1)T is a direction descent")
else
    disp("dk is not a direction descent")
end
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
ff=@(x,y) 5*x.^2+7*y.^2-3*x.*y;
dfdx=@(x,y) 10*x - 3*y;
dfdy=@(x,y) 14*y - 3*x;
%% Lambda calculation using Goldstein Armijo criteria
%% Xk+1 = Xk + lambda * dk
%% Choosing lambda > 0  
% Following criteria's are used to check the lambda > 0
% 1. f(Xk+1)-f(Xk) <= alpha*lambda * dk' * gradf(Xk)
% 2. dk'*gradf(Xk+1) >= eta * dk' * gradf(Xk)
%% Let initial lambda=0.1
%% Let initial guess x0 = [0;-1]
fcontour(ff,[-4 4 -2 10],'LevelStep',10)
grid
x0=0;y0=-1;
hold on
plot(x0,y0,'rx');
x0 =[0;-1];
d0=-[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))];
lam0=0.1;
x1=x0+lam0*d0;
alpha=0.0001;
eta=0.9;
i=0;
ITMAX=1000;
betaa1=100;
epsilon=0.00000001;
lam1=0.1;
tic
while (i < ITMAX) && (betaa1 >= epsilon)
    lam1=0.05+(0.05)*rand(1,1);
    if (ff(x1(1,1),x1(2,1))-ff(x0(1,1),x0(2,1))) <= (alpha*lam1*-[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))]'*[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))])
        if -[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))]'*[dfdx(x1(1,1),x1(2,1));dfdy(x1(1,1),x1(2,1))] >= (eta*-[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))]'*[dfdx(x0(1,1),x0(2,1));dfdy(x0(1,1),x0(2,1))])
            x2=x1+lam1*-[dfdx(x1(1,1),x1(2,1));dfdy(x1(1,1),x1(2,1))];
            plot(x2(1,1),x2(2,1),'rx');
            pause(0.05);
            betaa1=((dfdx(x2(1,1),x2(2,1))).^2+(dfdy(x2(1,1),x2(2,1))).^2)/(1+abs(ff(x2(1,1),x2(2,1))));
            x0=x1;
            x1=x2;
            i=i+1;
        end
    end
end
toc
hold off
disp("Local minimum solution approach: Cauchy's steepest descent method")
fprintf("The solution is obtained after %d iterations\n",i)
fprintf("The local minimum point : %0.5f\n",x1)
