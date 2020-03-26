%% Homework - 4
%% Submitted by Srinivasan Balan (200209659)
%% Problem 1
%% Apply the conjugate gradient method to solve the problem
%% State the problem
% Quadratic programming problem
% Idea to solve using Conjugate gradient as system of equations
%% State the solution
x =[1.0000;-2.0000;3.0000;0;0;1.0000]
%% State the method(s) used
% Used Conjugate gradient method
%% Show several starting and ending iterations (when applicable)
% Check the table after running code
%% Discuss any issues that occurred during implementation, if any
% Nothing. It run well
%% Present code
A=[4 0 0 1 0 0;0 4 0 0 1 0;0 0 5 0 0 1; 1 0 0 5 0 0; 0 1 0 0 6 0;0 0 1 0 0 6];
b=[4 ;-8; 16; 1 ;-2 ;9];
x=zeros(1,6);
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
tic
for i = 1:10
    if sqrt(r1(1,1)^2+r1(2,1)^2) < epsilon
        break
    elseif i > 0    
        beta1=r1'*r1*inv(rnot'*rnot);
        rnot=r1;
        p1=r1+beta1*pnot;
        pnot=p1;
        alpha1=r1'*r1 * inv(p1'*A*p1);
        x2=x1+alpha1*p1;
        x1=x2;
        r2=r1-alpha1*A*p1;
        r1=r2;
        iter=iter+1;
    end
end
toc
disp("Local minimum solution approach: Conjugate gradient method");
fprintf("The solution is obtained after %d iterations\n",iter);
fprintf("The values of x: %0.5f\n",x1);
%% Problem 2
%% Apply the non-linear conjugate gradient method to solve the problem
%% perform two iterations for n = 4 with x0=[1 1 1 1]'
%% State the problem
% As listed above, Non-linear equation to find the minimum 
% Idea to solve using Conjugate gradient as Fletcher Reeves method
%% State the solution
x = [1.5383;0.8770;.6622;0.5754]; % from Nelder Mead method
x = [0.86215;1.03801;0.88117;0.90059]; % terminated at 18 iterations
%% State the method(s) used
% Used Conjugate gradient Fletcher Reeves method
% Used linesearch as subroutine to find perfect line search for step length
% Coded Armijo search for additional linesearch method to find step length
%% Show several starting and ending iterations (when applicable)
% Check the table after running code
%% Discuss any issues that occurred during implementation, if any
% Majorly on step length convergence. After 18 iterations, step length
% never gave solution for F(Xk+1) < F(Xk)
%% Code
clear all
clc
syms x y z w
f=(x-2*y^2)^2+(y-2*z^2)^2+(z-2*w^2)^2;
gradient(f);
ff=@(x,y,z,w) (x-2*y.^2).^2+(y-2*z.^2).^2+(z-2*w.^2).^2;
dfdx=@(x,y,z,w) -8*w*(- 2*w^2 + z);
dfdy=@(x,y,z,w) - 4*y^2 + 2*x;
dfdz=@(x,y,z,w)  - 4*z^2 + 2*y - 8*y*(x - 2*y^2);
dfdw=@(x,y,z,w)  - 4*w^2 + 2*z - 8*z*(y - 2*z^2);
fff=@(x)(x(1)-2*x(2)^2)^2+(x(2)-2*x(3)^2)^2+(x(3)-2*x(4)^2)^2;
dfdxx=@(x) -8*x(4)*(- 2*x(4)^2 + x(3));
dfdyy=@(x) -4*x(2)^2 + 2*x(1);
dfdzz=@(x) - 4*x(3)^2 + 2*x(2) - 8*x(2)*(x(1) - 2*x(2)^2);
dfdww=@(x) - 4*x(3)^2 + 2*x(2) - 8*x(3)*(x(2) - 2*x(3)^2);
x=[1;1;1;1];
dfunc= [dfdxx(x);dfdyy(x);dfdzz(x);dfdww(x)];
xnot=x;
rnot=-[dfdx(x(1,1),x(2,1),x(3,1),x(4,1));dfdy(x(1,1),x(2,1),x(3,1),x(4,1));...
    dfdz(x(1,1),x(2,1),x(3,1),x(4,1));dfdw(x(1,1),x(2,1),x(3,1),x(4,1))];
pnot=rnot;
beta=0;
epsilon=0.00000001;
%% Choosing the step length lambda as low as possible lambda > 0
lam=linesearch(fff,dfunc,x); % Choose Initial step length of 0.01
%lam=goldarm(fff,dfdxx,dfdyy,dfdzz,dfdww,x); % Alternate step length
x1=xnot+lam*pnot;
iter=1;
betaa1=100;
fprintf("X1\t\t\tX2\t\t\tX3\t\t\tX4\t\t\tfun_Xk\n");
fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",xnot(1,1),xnot(2,1),xnot(3,1),xnot(4,1),ff(xnot(1,1),xnot(2,1),xnot(3,1),xnot(4,1)));
fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x1(1,1),x1(2,1),x1(3,1),x1(4,1),ff(x1(1,1),x1(2,1),x1(3,1),x1(4,1)));
tic
for i = 1:18
    % I have provided the maximum runs of iteration till step length
    % explodes and f(x) starts increasing
    % We can keep "while ITMAX < 3" as a condition to keep 2 iterations
    % by initializing ITMAX=0 to begin with. 
    if (betaa1 < epsilon) 
        break
    elseif i > 0    
        beta1=[dfdx(x1(1,1),x1(2,1),x1(3,1),x1(4,1));dfdy(x1(1,1),x1(2,1),x1(3,1),x1(4,1));...
    dfdz(x1(1,1),x1(2,1),x1(3,1),x1(4,1));dfdw(x1(1,1),x1(2,1),x1(3,1),x1(4,1))]'*[dfdx(x1(1,1),x1(2,1),x1(3,1),x1(4,1));dfdy(x1(1,1),x1(2,1),x1(3,1),x1(4,1));...
    dfdz(x1(1,1),x1(2,1),x1(3,1),x1(4,1));dfdw(x1(1,1),x1(2,1),x1(3,1),x1(4,1))]*inv([dfdx(x(1,1),x(2,1),x(3,1),x(4,1));dfdy(x(1,1),x(2,1),x(3,1),x(4,1));...
    dfdz(x(1,1),x(2,1),x(3,1),x(4,1));dfdw(x(1,1),x(2,1),x(3,1),x(4,1))]'*[dfdx(x(1,1),x(2,1),x(3,1),x(4,1));dfdy(x(1,1),x(2,1),x(3,1),x(4,1));...
    dfdz(x(1,1),x(2,1),x(3,1),x(4,1));dfdw(x(1,1),x(2,1),x(3,1),x(4,1))]);
        p1=-[dfdx(x(1,1),x(2,1),x(3,1),x(4,1));dfdy(x(1,1),x(2,1),x(3,1),x(4,1));...
    dfdz(x(1,1),x(2,1),x(3,1),x(4,1));dfdw(x(1,1),x(2,1),x(3,1),x(4,1))]+beta1*pnot;
        pnot=p1;
        dfunc=[dfdxx(x1);dfdyy(x1);dfdzz(x1);dfdww(x1)];
        lam1=linesearch(fff,dfunc,x1);
        %lam=goldarm(fff,dfdxx,dfdyy,dfdzz,dfdww,x1); % Alternate step length
        x2=x1+lam1*p1;
        betaa1=dfdx(x2(1,1),x2(2,1),x2(3,1),x2(4,1)).^2+dfdy(x2(1,1),x2(2,1),x2(3,1),x2(4,1)).^2+dfdz(x2(1,1),x2(2,1),x2(3,1),x2(4,1)).^2+dfdw(x2(1,1),x2(2,1),x2(3,1),x2(4,1)).^2;
        x=x1;
        x1=x2;
        iter=iter+1;
        fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x2(1,1),x2(2,1),x2(3,1),x2(4,1),ff(x2(1,1),x2(2,1),x2(3,1),x2(4,1)));
    end
end
toc
disp("Local minimum solution approach: Non-Linear Conjugate gradient method");
fprintf("Solution is terminated after %d iterations with step length 0.0001\n",iter+1);
fprintf("The values of x: %0.5f\n",x1);
%% Using fminsearch to get the solution when n=4 with starting solution =[1 1 1 1]
f=@(x)(x(1)-2*x(2)^2)^2+(x(2)-2*x(3)^2)^2+(x(3)-2*x(4)^2)^2
[x,feval]=fminsearch(f,[1 1 1 1])
% x = [1.5383;0.8770;.6622;0.5754] as solution from Nelder Mead method
