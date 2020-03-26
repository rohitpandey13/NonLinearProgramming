%% Homework - 6
%% Submitted by Srinivasan Balan (200209659)
%% Solve the problem to minimize using the following methods
clear all
clc
ff=@(x) x(1)^2+x(2)^2-6*x(1)-8*x(2)+10;
con1=@(x) 4*x(1)^2+x(2)^2-16;
con2=@(x) 3*x(1)+5*x(2)-4;
%con3=@(x) -x(1);
%con4=@(x) -x(2);
%% a. Exterior Penalty method





%% b. Interior penalty method





%% c. Zoutendijk Method







%% Using MATLAB fmincon function
clear all
clc
format compact
ff=@(x) x(1)^2+x(2)^2-6*x(1)-8*x(2)+10;
x0=[0;0];
lb=[0;0];
ub=[Inf; Inf];
x=fmincon(ff,x0,[],[],[],[],lb,ub,@nlcon);
fprintf('Initial objective: %d\n',ff(x0))
fprintf('The values of x: %0.5f\n',x)
fprintf('Final objective: %d\n',ff(x))

 
