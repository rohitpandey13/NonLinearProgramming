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
%dfunc= [dfdxx(x);dfdyy(x);dfdzz(x);dfdww(x)];
xnot=x;
rnot=-[dfdx(x(1,1),x(2,1),x(3,1),x(4,1));dfdy(x(1,1),x(2,1),x(3,1),x(4,1));...
    dfdz(x(1,1),x(2,1),x(3,1),x(4,1));dfdw(x(1,1),x(2,1),x(3,1),x(4,1))];
pnot=rnot;
beta=0;
epsilon=0.00000001;
%% Choosing the step length lambda as low as possible lambda > 0
lam=goldarm(fff,dfdxx,dfdyy,dfdzz,dfdww,x); % Choose Initial step length of 0.01
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
        %dfunc=[dfdxx(x1);dfdyy(x1);dfdzz(x1);dfdww(x1)];
        lam1=goldarm(fff,dfdxx,dfdyy,dfdzz,dfdww,x1);
        x2=x1+lam1*p1;
        betaa1=dfdx(x2(1,1),x2(2,1),x2(3,1),x2(4,1)).^2+dfdy(x2(1,1),x2(2,1),x2(3,1),x2(4,1)).^2+dfdz(x2(1,1),x2(2,1),x2(3,1),x2(4,1)).^2+dfdw(x2(1,1),x2(2,1),x2(3,1),x2(4,1)).^2;
        x=x1;
        x1=x2;
        iter=iter+1;
        fprintf("%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.5f\n",x2(1,1),x2(2,1),x2(3,1),x2(4,1),ff(x2(1,1),x2(2,1),x2(3,1),x2(4,1)));
    end
end
toc
disp("Local minimum solution approach: Non-Linear Conjugate gradient method")
fprintf("Solution is terminated after %d iterations with step length 0.0001\n",iter+1)
fprintf("The values of x: %0.5f\n",x1)