%% Homework - 5
%% Submitted by Srinivasan Balan (200209659)
%% Problem 1
clear all
clc
syms x y z
f=0.5*(x^2+5*y^2+25*z^2)+x+y+z
gradient(f)
A=[1 0 0;0 5 0;0 0 25]
b=[-1;-1;-1]
x=zeros(1,3);
ff=@(x,y,z) 0.5*(x.^2+5*y.^2+25*z.^2)+x+y+z;
dfdx=@(x,y,z) x+1;
dfdy=@(x,y,z) 5*y + 1;
dfdz=@(x,y,z) 25*z + 1;
%% Apply the conjugate gradient method to solve the problem
xnot=x';
% Minimizing the quadratic form is equivalent to solving Ax=b
% system of linear equations
% Faster than steepest descent and takes only n iterations as number of
% dimensions of the A matrix. This problem converges in 6 iterations. 
Ex=inv(A)*b % Solution to Ax=b
% x=[1;-2;3;0;0;1] as solution to this problem
rnot=b-A*xnot;
pnot=rnot;
epsilon=0.00000001;
alphanot=rnot'*rnot * inv(pnot'*A*pnot);
x1=xnot+alphanot*pnot;
r1=rnot-alphanot*A*pnot;
betanot=r1'*r1*inv(rnot'*rnot);
p1=r1+betanot*pnot;
iter=1;
betaa1=100;
fprintf("X1\t\t\tX2\t\t\tX3\t\t\tfun_Xk\n");
fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",xnot(1,1),xnot(2,1),xnot(3,1),ff(xnot(1,1),xnot(2,1),xnot(3,1)));
fprintf("%0.5f\t%0.5f\t%0.5f\t%0.5f\n",x1(1,1),x1(2,1),x1(3,1),ff(x1(1,1),x1(2,1),x1(3,1)));
tic
for i = 1:10
    if (sqrt(r1(1,1)^2+r1(2,1)^2) < epsilon) || (betaa1 < epsilon) 
        break
    elseif i > 0    
        beta1=r1'*r1*inv(rnot'*rnot);
        rnot=r1;
        p1=r1+beta1*pnot;
        pnot=p1;
        alpha1=r1'*r1 * inv(p1'*A*p1);
        x2=x1+alpha1*p1;
        betaa1=dfdx(x2(1,1),x2(2,1),x2(3,1)).^2+dfdy(x2(1,1),x2(2,1),x2(3,1)).^2+dfdz(x2(1,1),x2(2,1),x2(3,1)).^2;
        x1=x2;
        r2=r1-alpha1*A*p1;
        r1=r2;
        iter=iter+1;
        fprintf("%0.5f\t%0.5f\t%0.5f\t%0.5f\n",x2(1,1),x2(2,1),x2(3,1),ff(x2(1,1),x2(2,1),x2(3,1)));
    end
end
toc
disp("Local minimum solution approach: Conjugate gradient method");
fprintf("The solution is obtained after %d iterations\n",iter);
fprintf("The values of x: %0.5f\n",x1);
fprintf("The objective function value: %0.5f\n",ff(x1(1,1),x1(2,1),x1(3,1)));
%% Problem 2
%% Using BFGS method in Quasi Newton approach to find the minimum of a function
%% Using I as initial Hessian approximation
clear all
clc
ITEMAX=1000;
epsilon=0.00000001;
%% 1. Problem
%syms x y
%f=x^2+2*y^2-2*x*y-2*y;
%gradient(f);
f1=@(x) x(1)^2+2*x(2)^2-2*x(1)*x(2)-2*x(2);
df1x=@(x) 2*x(1) - 2*x(2);
df1y=@(x) 4*x(2) - 2*x(1) - 2;
xnot=[0;0];
dfunc=[df1x(xnot);df1y(xnot)];
%% BFGS code
I=eye(2,2);
H0=I;
d0=-H0*dfunc;
lam=linesearch(f1,dfunc,xnot); % Choose Initial step length of 0.01
%lam=goldarm(f1,df1x,df1y,xnot); % Alternate step length
x1=xnot+lam*d0;
iter=1;
betaa1=100;
fprintf("X1\t\t\tX2\t\t\tStepLength\tdirectionX1\tdirectionX2\tfun_Xk\n");
fprintf("%0.5f\t\t%0.5f\t\t\t\t\t\t\t\t\t\t\t%0.5f\n",xnot(1,1),xnot(2,1),f1(xnot));
fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x1(1),x1(2),lam,d0(1),d0(2),f1(x1));
tic
S0=lam*d0;
y0=[df1x(x1);df1y(x1)]-[df1x(xnot);df1y(xnot)];
H1=H0+((S0-H0*y0)*S0' + S0*(S0-H0*y0)')*inv(y0'*S0) - (((S0-H0*y0)'*y0)*(S0*S0'))*(inv(y0'*S0))^2;                
tic
for i = 1:10000
    % I have provided the maximum runs of iteration till step length
    % explodes and f(x) starts increasing
    % Function line search is used for step length calculation
    % We can keep "while ITMAX < 1000" as a condition
    % by initializing iter=1 to begin with as mentioned below
    if (betaa1 < epsilon) || (iter > ITEMAX)
        break
    elseif i > 0    
        d1=-H1*[df1x(x1);df1y(x1)];
        dfunc=[df1x(x1);df1y(x1)];
        lam1=linesearch(f1,dfunc,x1);
        %lam1=goldarmb(f1,df1x,df1y,x1); % Alternate step length
        x2=x1+lam1*d1;
        fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x2(1),x2(2),lam1,d1(1),d1(2),f1(x2));
        betaa1=((df1x(x2)).^2+(df1y(x2)).^2)/(1+abs(f1(x2)));
        S1=lam1*d1;
        y1=[df1x(x2);df1y(x2)]-[df1x(x1);df1y(x1)];
        x1=x2;
        H2=H1 + ((S1-H1*y1)*S1' + S1*(S1-H1*y1)')*inv(y1'*S1) - (((S1-H1*y1)'*y1)*(S1*S1'))*(inv(y1'*S1))^2;
        H1=H2;
        iter=iter+1;
    end
end
toc
disp("Local minimum solution approach: BFGS in Quasi Newton method");
fprintf("Solution is terminated after %d iterations with step length 0.0001 at each step\n",iter+1);
fprintf("The values of x: %0.5f\n",x1);
fprintf("The value of the objective: %0.5f\n",f1(x1));
%% Solution using Nelder Mead search
[x,feval]=fminsearch(f1,[0 0])
%% 2. Problem
%% c=1
%% Using I as initial Hessian approximation
clear all
clc
ITEMAX=1000;
epsilon=0.00000001;
%syms x y
%fc=(x-1)^2+(y-1)^2+(x^2+y^2-0.25)^2;
%gradient(fc)
f1=@(x) (x(1)-1)^2+(x(2)-1)^2+(x(1)^2+x(2)^2-0.25)^2;
df1x=@(x) 2*x(1) + 4*x(1)*(x(1)^2 + x(2)^2 - 1/4) - 2;
df1y=@(x)  2*x(2) + 4*x(2)*(x(1)^2 + x(2)^2 - 1/4) - 2;
xnot=[1;-1];
dfunc=[df1x(xnot);df1y(xnot)];
%% BFGS code
I=eye(2,2);
H0=I;
d0=-H0*dfunc;
lam=linesearch(f1,dfunc,xnot); % Choose Initial step length of 0.01
%lam=goldarm(f1,df1x,df1y,xnot); % Alternate step length
x1=xnot+lam*d0;
iter=1;
betaa1=100;
fprintf("X1\t\t\tX2\t\t\tStepLength\tdirectionX1\tdirectionX2\tCond\tfun_Xk\n");
fprintf("-------------------------------------------------------------------------\n");
fprintf("%0.5f\t\t%0.5f\t\t\t\t\t\t\t\t\t\t\t\t%0.5f\n",xnot(1,1),xnot(2,1),f1(xnot));
fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x1(1),x1(2),lam,d0(1),d0(2),cond(H0),f1(x1));
tic
S0=lam*d0;
y0=[df1x(x1);df1y(x1)]-[df1x(xnot);df1y(xnot)];
H1=H0+((S0-H0*y0)*S0' + S0*(S0-H0*y0)')*inv(y0'*S0) - (((S0-H0*y0)'*y0)*(S0*S0'))*(inv(y0'*S0))^2;                
tic
for i = 1:10000
    % I have provided the maximum runs of iteration till step length
    % explodes and f(x) starts increasing
    % We can keep "while ITMAX < 1000" as a condition
    % by initializing iter=1 to begin with as mentioned below
    if (betaa1 < epsilon) || (iter > ITEMAX)
        break
    elseif i > 0    
        d1=-H1*[df1x(x1);df1y(x1)];
        dfunc=[df1x(x1);df1y(x1)];
        lam1=linesearch(f1,dfunc,x1);
        %lam1=goldarmb(f1,df1x,df1y,x1); % Alternate step length
        x2=x1+lam1*d1;
        fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x2(1),x2(2),lam1,d1(1),d1(2),cond(H1),f1(x2));
        betaa1=((df1x(x2)).^2+(df1y(x2)).^2)/(1+abs(f1(x2)));
        S1=lam1*d1;
        y1=[df1x(x2);df1y(x2)]-[df1x(x1);df1y(x1)];
        x1=x2;
        H2=H1 + ((S1-H1*y1)*S1' + S1*(S1-H1*y1)')*inv(y1'*S1) - (((S1-H1*y1)'*y1)*(S1*S1'))*(inv(y1'*S1))^2;
        H1=H2;
        iter=iter+1;
    end
end
toc
disp("Local minimum solution approach: BFGS in Quasi Newton method");
fprintf("Solution is terminated after %d iterations with step length 0.0001 at each step\n",iter+1);
fprintf("The values of x: %0.5f\n",x1);
fprintf("The value of the objective: %0.5f\n",f1(x1));
%% Solution using Nelder Mead search
[x,feval]=fminsearch(f1,[1 -1])
%% c=10
clear all
clc
ITEMAX=1000;
epsilon=0.00000001;
syms x y
fc2=(x-1)^2+(y-1)^2+10*(x^2+y^2-0.25)^2;
gradient(fc2);
f1=@(x) (x(1)-1)^2+(x(2)-1)^2+10*(x(1)^2+x(2)^2-0.25)^2;
df1x=@(x) 2*x(1) + 40*x(1)*(x(1)^2 + x(2)^2 - 1/4) - 2;
df1y=@(x) 2*x(2) + 40*x(2)*(x(1)^2 + x(2)^2 - 1/4) - 2;
xnot=[1;-1];
dfunc=[df1x(xnot);df1y(xnot)];
%% BFGS code
I=eye(2,2);
H0=I;
d0=-H0*dfunc;
lam=linesearch(f1,dfunc,xnot); % Choose Initial step length of 0.01
%lam=goldarm(f1,df1x,df1y,xnot); % Alternate step length
x1=xnot+lam*d0;
iter=1;
betaa1=100;
fprintf("X1\t\t\tX2\t\t\tStepLength\tdirectionX1\tdirectionX2\tCond\tfun_Xk\n");
fprintf("-------------------------------------------------------------------------\n");
fprintf("%0.5f\t\t%0.5f\t\t\t\t\t\t\t\t\t\t\t\t%0.5f\n",xnot(1,1),xnot(2,1),f1(xnot));
fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x1(1),x1(2),lam,d0(1),d0(2),cond(H0),f1(x1));
tic
S0=lam*d0;
y0=[df1x(x1);df1y(x1)]-[df1x(xnot);df1y(xnot)];
H1=H0+((S0-H0*y0)*S0' + S0*(S0-H0*y0)')*inv(y0'*S0) - (((S0-H0*y0)'*y0)*(S0*S0'))*(inv(y0'*S0))^2;                
tic
for i = 1:10000
    % I have provided the maximum runs of iteration till step length
    % explodes and f(x) starts increasing
    % We can keep "while ITMAX < 1000" as a condition
    % by initializing iter=1 to begin with as mentioned below
    if (betaa1 < epsilon) || (iter > ITEMAX)
        break
    elseif i > 0    
        d1=-H1*[df1x(x1);df1y(x1)];
        dfunc=[df1x(x1);df1y(x1)];
        lam1=linesearch(f1,dfunc,x1);
        %lam1=goldarmb(f1,df1x,df1y,x1); % Alternate step length
        x2=x1+lam1*d1;
        fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x2(1),x2(2),lam1,d1(1),d1(2),cond(H1),f1(x2));
        betaa1=((df1x(x2)).^2+(df1y(x2)).^2)/(1+abs(f1(x2)));
        S1=lam1*d1;
        y1=[df1x(x2);df1y(x2)]-[df1x(x1);df1y(x1)];
        x1=x2;
        H2=H1 + ((S1-H1*y1)*S1' + S1*(S1-H1*y1)')*inv(y1'*S1) - (((S1-H1*y1)'*y1)*(S1*S1'))*(inv(y1'*S1))^2;
        H1=H2;
        iter=iter+1;
    end
end
toc
disp("Local minimum solution approach: BFGS in Quasi Newton method");
fprintf("Solution is terminated after %d iterations with step length 0.0001 at each step\n",iter+1);
fprintf("The values of x: %0.5f\n",x1);
fprintf("The value of the objective: %0.5f\n",f1(x1));
%% Solution using Nelder Mead search
[x,feval]=fminsearch(f1,[1 -1])
%% c=100
clear all
clc
ITEMAX=1000;
epsilon=0.00000001;
syms x y
fc3=(x-1)^2+(y-1)^2+100*(x^2+y^2-0.25)^2;
gradient(fc3);
f1=@(x) (x(1)-1)^2+(x(2)-1)^2+100*(x(1)^2+x(2)^2-0.25)^2;
df1x=@(x) 2*x(1) + 400*x(1)*(x(1)^2 + x(2)^2 - 1/4) - 2;
df1y=@(x) 2*x(2) + 400*x(2)*(x(1)^2 + x(2)^2 - 1/4) - 2;
xnot=[1;-1];
dfunc=[df1x(xnot);df1y(xnot)];
%% BFGS code
I=eye(2,2);
H0=I;
d0=-H0*dfunc;
lam=linesearch(f1,dfunc,xnot); % Choose Initial step length of 0.01
%lam=goldarm(f1,df1x,df1y,xnot); % Alternate step length
x1=xnot+lam*d0;
iter=1;
betaa1=100;
fprintf("X1\t\t\tX2\t\t\tStepLength\tdirectionX1\tdirectionX2\tCond\tfun_Xk\n");
fprintf("-------------------------------------------------------------------------\n");
fprintf("%0.5f\t\t%0.5f\t\t\t\t\t\t\t\t\t\t\t\t%0.5f\n",xnot(1,1),xnot(2,1),f1(xnot));
fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x1(1),x1(2),lam,d0(1),d0(2),cond(H0),f1(x1));
tic
S0=lam*d0;
y0=[df1x(x1);df1y(x1)]-[df1x(xnot);df1y(xnot)];
H1=H0+((S0-H0*y0)*S0' + S0*(S0-H0*y0)')*inv(y0'*S0) - (((S0-H0*y0)'*y0)*(S0*S0'))*(inv(y0'*S0))^2;                
tic
for i = 1:10000
    % I have provided the maximum runs of iteration till step length
    % explodes and f(x) starts increasing
    % We can keep "while ITMAX < 1000" as a condition
    % by initializing iter=1 to begin with as mentioned below
    if (betaa1 < epsilon) || (iter > ITEMAX)
        break
    elseif i > 0    
        d1=-H1*[df1x(x1);df1y(x1)];
        dfunc=[df1x(x1);df1y(x1)];
        lam1=linesearch(f1,dfunc,x1);
        %lam1=goldarmb(f1,df1x,df1y,x1); % Alternate step length
        x2=x1+lam1*d1;
        fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x2(1),x2(2),lam1,d1(1),d1(2),cond(H1),f1(x2));
        betaa1=((df1x(x2)).^2+(df1y(x2)).^2)/(1+abs(f1(x2)));
        S1=lam1*d1;
        y1=[df1x(x2);df1y(x2)]-[df1x(x1);df1y(x1)];
        x1=x2;
        H2=H1 + ((S1-H1*y1)*S1' + S1*(S1-H1*y1)')*inv(y1'*S1) - (((S1-H1*y1)'*y1)*(S1*S1'))*(inv(y1'*S1))^2;
        H1=H2;
        iter=iter+1;
    end
end
toc
disp("Local minimum solution approach: BFGS in Quasi Newton method");
fprintf("Solution is terminated after %d iterations with step length 0.0001 at each step\n",iter+1);
fprintf("The values of x: %0.5f\n",x1);
fprintf("The value of the objective: %0.5f\n",f1(x1));
%% Solution using Nelder Mead search
[x,feval]=fminsearch(f1,[1 -1])
%% Effect on algorithm performance using Condition number
% Condition number estimates worst case loss of precision
% How much the output value of the function can change for a 
% small change in the input argument.
% When C=1,10, the problem is well conditioned with low condition number
% when C=100, the problem becomes ill-conditioned and it increases the
% number of iterations to more than 1000. The run time is around 3 seconds.
%% Problem 3
%% Nelder Mead method in a closed ball of radius 10
%% 1. Problem
clear all
clc
ff=@(x) x(1)-x(2)+2*x(1)^2+2*x(1)*x(2)
r=10;
%% Selecting the initial simplex
rand('state',123);
n=3;
done=false;
while ~done
XY=rand(n,2);
XY=XY*2;
isin=sum(XY.^2,2)<r^2;
m=sum(isin);
if m==3
    done=true;
    break
end
end
XY
%% Initialization of parameters
%% ANMS Implementation of Nelder Mead algorithm
% Defined inside the function NMead for reference.
% As per ANMS implementation: The following parameters give better results
% n represents the number of dimensions of the problem
% rho=1;
% xi=1+2/n;
% gam=0.75-1/2n;
% sig=1-1/n;
%% Algorithm implementation with hyper parameters
% Remove the comment given in the 'pause' for contour plotting at 
% line 40 and 135 in the NMead function
% max_feval = 1000; Preset the number of ITEMAX.
% If ITEMAX is more than 1000, still convergence not reached, it is assumed
% to be diverged. 
% This will be printed with some Nelder Mead - Warning signal !
% Flag is used to plot contour for n=2 dimensions problem only.
%% Solution using NMead function
[x_opt,n_feval]=NMead(XY,ff,1);
ff(x_opt);
disp("Local minimum solution approach: ANMS Implementation using Nelder Mead method");
fprintf("Solution is obtained after %d iterations\n",n_feval);
fprintf("The values of x: %0.5f\n",x_opt);
fprintf("The value of the objective: %0.5f\n",ff(x_opt));
disp("This problem is not converging to the minimum point");
%% Matlab Nelder Mead method solution
[x,feval]=fminsearch(ff,XY(1,:))
%% 2. Problem
clear all
clc
fg=@(x) x(1)^2+2*x(2)^2-2*x(1)*x(2)-2*x(2)
r=10;
%% Selecting the initial simplex
rand('state',123);
n=3;
done=false;
while ~done
XY=rand(n,2);
XY=XY*2;
isin=sum(XY.^2,2)<r^2;
m=sum(isin);
if m==3
    done=true;
    break
end
end
XY
%% Initialization of parameters
%% ANMS Implementation of Nelder Mead algorithm
% Defined inside the function NMead for reference.
% As per ANMS implementation: The following parameters give better results
% n represents the number of dimensions of the problem
% rho=1;
% xi=1+2/n;
% gam=0.75-1/2n;
% sig=1-1/n;
% Defined inside the function NMead for reference. 
%% Algorithm implementation with hyper parameters
% Remove the comment given in the 'pause' for contour plotting at 
% line 40 and 135 in the NMead function
% max_feval = 1000; Preset the number of ITEMAX.
% If ITEMAX is more than 1000, still convergence not reached, it is assumed
% to be diverged.
% This will be printed with some Nelder Mead - Warning signal !
% Flag is used to plot contour for n=2 dimensions problem only. 
%% Solution using NMead function
[x_opt,n_feval]=NMead(XY,fg,1);
fg(x_opt);
disp("Local minimum solution approach: ANMS Implementation using Nelder Mead method");
fprintf("Solution is obtained after %d iterations\n",n_feval);
fprintf("The values of x: %0.5f\n",x_opt);
fprintf("The value of the objective: %0.5f\n",fg(x_opt));
%% Matlab Nelder Mead method solution
[x,feval]=fminsearch(fg,XY(1,:))


