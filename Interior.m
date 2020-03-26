%% Homework - 6
%% Submitted by Srinivasan Balan (200209659)
%% Solve the problem to minimize using the following methods
clear all
clc
con1=@(x) 4*x(1)^2+x(2)^2-16;
con2=@(x) 3*x(1)+5*x(2)-4;
%con3=@(x) -x(1);
%con4=@(x) -x(2);
%% a. Interior Penalty method
ITEMAX=1000;
epsilon=0.00000001;
%nmin=0.000000001;
k=0;
ite=0;
tic
for mu=[10^1 10^0 10^-1 10^-2 10^-3 10^-4 10^-5 10^-6 10^-7 10^-8 10^-9 10^-10]
%% Log function
f1=@(x) x(1)^2+x(2)^2-6*x(1)-8*x(2)+10-mu*log((4*x(1)^2+x(2)^2-16))-mu*log((3*x(1)+5*x(2)-4))-mu*log(x(1)^2)-mu*log(x(2)^2);
df1x=@(x) 2*x(1) - (3*mu)/(3*x(1) + 5*x(2) - 4) - (2*mu)/x(1) - (8*x(1)*mu)/(4*x(1)^2 + x(2)^2 - 16) - 6;
df1y=@(x) 2*x(2) - (5*mu)/(3*x(1) + 5*x(2) - 4) - (2*mu)/x(2) - (2*x(2)*mu)/(4*x(1)^2 + x(2)^2 - 16) - 8;
uu=1;
xx=[0.5;0.5];
while uu>0
if con1(xx)<0 && con2(xx)<0
    uu=0;
    xnot=xx;
else
    xx=[-rand; -rand];
end
end
dfunc=[df1x(xnot);df1y(xnot)];

%% Inverse function
%f1=@(x) x(1)^2+x(2)^2-6*x(1)-8*x(2)+10-mu*((4*x(1)^2+x(2)^2-16)^-1)-mu*((3*x(1)+5*x(2)-4)^-1)-mu*(x(1)^-2)-mu*(x(2)^-2);
%df1x=@(x)  2*x(1) + (3*mu)/(3*x(1) + 5*x(2) - 4)^2 + (2*mu)/x(1)^3 + (8*x(1)*mu)/(4*x(1)^2 + x(2)^2 - 16)^2 - 6;
%df1y=@(x)  2*x(2) + (5*mu)/(3*x(1) + 5*x(2) - 4)^2 + (2*mu)/x(2)^3 + (2*x(2)*mu)/(4*x(1)^2 + x(2)^2 - 16)^2 - 8;
%uu=1;
%xx=[-1;-1];
%while uu>0
%if con1(xx)<0 && con2(xx)<0
 %   uu=0;
  %  xnot=xx;
%else
 %   xx=[-rand; -rand];
%end
%end
%dfunc=[df1x(xnot);df1y(xnot)];
%% BFGS code
I=eye(2,2);
H0=I;
d0=-H0*dfunc;
lam=linesearch1(f1,dfunc,xnot); % Choose Initial step length of 0.01
%lam=goldarm(f1,df1x,df1y,xnot); % Alternate step length
x1=xnot+lam*d0;
iter=1;
betaa1=100;
S0=lam*d0;
y0=[df1x(x1);df1y(x1)]-[df1x(xnot);df1y(xnot)];
H1=H0+((S0-H0*y0)*S0' + S0*(S0-H0*y0)')*inv(y0'*S0) - (((S0-H0*y0)'*y0)*(S0*S0'))*(inv(y0'*S0))^2;                
for i = 1:10000
    if (betaa1 < epsilon) || (iter > ITEMAX)
        break
    elseif i > 0    
        d1=-H1*[df1x(x1);df1y(x1)];
        dfunc=[df1x(x1);df1y(x1)];
        lam1=linesearch1(f1,dfunc,x1);
        x2=x1+lam1*d1;
        betaa1=((df1x(x2)).^2+(df1y(x2)).^2)/(1+abs(f1(x2)));
        S1=lam1*d1;
        y1=[df1x(x2);df1y(x2)]-[df1x(x1);df1y(x1)];
        x1=x2;
        H2=H1 + ((S1-H1*y1)*S1' + S1*(S1-H1*y1)')*inv(y1'*S1) - (((S1-H1*y1)'*y1)*(S1*S1'))*(inv(y1'*S1))^2;
        H1=H2;
        iter=iter+1;
    end
end
fprintf("Iteration : %d\n",ite);
fprintf("The values of x: %0.5f\n",x1);
fprintf("The value of the objective: %0.5f\n",f1(x1));
fprintf("_____________________________________________\n");
mu;
ite=ite+1;
end
toc
disp("Local minimum solution approach: BFGS in Quasi Newton method");
disp("Local minimum solution for Constrained problem: Interior penalty method");
fprintf("Solution is terminated after %d iterations\n",ite);
fprintf("The values of x: %0.5f\n",x1);
fprintf("The value of the objective: %0.5f\n",f1(x1));