%% OR 506 HW_2
%% 1. Find the minimum value of f(x)=e^(-x) + x^2
%% Optimal solution
clear all
clc
p = @(x) exp(-x)+x.*x
x=-100:0.1:100;
plot(x,p(x))
grid on
hold on
Z=fminsearch(p,1.1)
plot(Z,p(Z),'r*')
hold off
%% Golden section method
clc
clear all
func=@(x)exp(-x)+x.*x;
low=-10;
high=10;
%fplot(func,[low,high])
tol=0.00000001;
% Algorithm 1
goldsection(func,low,high,tol); 
% Algorithm 2
gold_section(func,low,high,tol);
%% Bisection method
clc
clear all
funct=@(x)exp(-x)+x.*x;
low=-100;
high=200;
tol=0.00000001;
bisection(funct,low,high,tol);
%% Fibonacci search method
clc
clear all
fun=@(x)exp(-x)+x.*x;
low=-10;
high=10;
tol=0.00000001;
Fibonacci(fun,low,high,tol);
%% 2. Newton's method to minimize f(x)=ln(e^x+e^-x)
%% Plotting the given function and solving using fminsearch
p=@(x) log(exp(-x)+exp(x))
x=-100:0.1:100;
plot(x,p(x))
grid on
hold on
Z=fminsearch(p,1.1)
plot(Z,p(Z),'r*')
hold off
%% Newton's method with starting point x0=1.1
p=@(x) log(exp(-x)+exp(x))
low=-100;
high=100;
plow=p(low)
phigh=p(high)
syms f(x)
f(x)=p;
deriv=diff(f);
delta=1000;
epsilon = 0.0001;
iter=0;
x=1.1; % Starting point x0=1.1; try different points close to 0
tic
while (delta > epsilon)
        y=x-(p(x)/double(subs(deriv(x),x)))
        delta=abs(x-y);
        x=y;
        iter=iter+1;
end
xopt=x;toc
disp("Local minimum solution approach: Newton's method")
fprintf("The solution is obtained after %d iterations\n",iter)
fprintf("The local minimum point : %0.5f\n",xopt)
% This function is not converging to 0 as it is asymptotic towards 
% infinity after crossing 0. The gradient will push or oscilate 
% between positive and negative numbers as the tangent will not reach 0
% Tried with different starting points and all points converging
% at infinity as Nan
% Therefore Newtons method will cycle back and forth between two
% value and never converge at all. The answer is at x=0.0
% This shows the failure of Newtons method on convergence
%% 3. Newton's method for solving equations
%% a. f(x) = 7x4+3x3+2x2+9x+4=0
%% Optimal solution using Syms
syms x
eqn = 7*x^4+3*x^3+2*x^2+9*x+4==0
sol=solve(eqn)
% Gives output in syms
%% Alitar (4 roots of the equation)
p=[7 3 2 9 4];
r=roots(p)
% Gives output in complex number format
%% Plotting function and solution using fzero
p=@(x) 7*x.^4+3*x.^3+2*x.^2+9*x+4
x=-100:0.1:100;
plot(x,p(x))
grid on
hold on
Z=fzero(p,0.5)
plot(Z,p(Z),'r*')
hold off
%% Newton's method
p=@(x) 7*x.^4+3*x.^3+2*x.^2+9*x+4;
epsilon=0.00000001;
syms f(x)
f(x)=p;
dev=diff(f);
iter=0;
x=1;
delta=1000;
tic
while (delta > epsilon)
     y=x-(p(x)/eval(dev(x)));
     delta=abs(y-x);
     x=y;
     iter=iter+1;
end
xopt=x;
toc
disp("Local minimum solution approach: Newton's method")
fprintf("The solution is obtained after %d iterations\n",iter)
fprintf("The roots of the equation : %0.5f\n",xopt)
%% b. 
%% f(x1,x2) = 3x1x2+7x1+2x2-3=0
%% f(x1,x2) = 5x1x2-9x1-4x2+6=0
%% Solving the non-linear system of equations
%% refer the roots.m file for the two functions
fun=@roots;
x0=[1,1];
options=optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
x=fsolve(fun,x0,options)
%% Using Fsolve method
F=@(x) [3*x(1)*x(2)+7*x(1)+2*x(2)-3;5*x(1)*x(2)-9*x(1)-4*x(2)+6];
Initial=[-1;1];
Options=optimset('Display','iter');
XY=fsolve(F,Initial,Options);
XY
sbzero=F(XY)
ezplot('3*x1*x2+7*x1+2*x2-3');
hold on
ezplot('5*x1*x2-9*x1-4*x2+6');
plot(XY(1),XY(2),'ro')
%% Newton's method
p1=@(x) 3*x(1)*x(2)+7*x(1)+2*x(2)-3;
p2=@(x) 5*x(1)*x(2)-9*x(1)-4*x(2)+6;
syms x1 x2;
F=[3*x1*x2+7*x1+2*x2-3;5*x1*x2-9*x1-4*x2+6];
x1=-1;
x2=1;
eval(F);
eval(jacobian(F));
iter=0;
x=[x1;x2];
u=100;
epsilon=0.001;
tic
while (u > epsilon)
     x1=x(1,1);
     x2=x(2,1);
     a=inv(eval(jacobian(F)));
     b=eval(F);
     y=x-a*b;
     delta=abs(y-x);
     u=norm(delta);
     x=y;
     iter=iter+1;
end
xopt=x;
toc
disp("Local minimum solution approach: Newton's method")
fprintf("The solution is obtained after %d iterations\n",iter)
fprintf("The roots of the equation : %0.5f\n",xopt) 












