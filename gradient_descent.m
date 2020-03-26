%% Use the gradient to x and y from gradient function
%% Steepest descent
clear all
clc
b=@(x,y)5*x.^2+7*y.^2-2*x.*y
dfdx=@(x,y) 10*x - 2*y
dfdy=@(x,y) 14*y - 2*x
%% Subplot shows contour
%subplot 211
fcontour(b,[-4 4 -2 10],'LevelStep',10)
grid
x0=0;y0=-2;
hold on
plot(x0,y0,'o');
i=0;
iter=0;
ITMAX=10000;
epsilon=0.00000001;
beta1=100;
beta2=100;tic
while i < ITMAX && (beta1 > epsilon || beta2 > epsilon)
    s1=-dfdx(x0,y0);
    s2=-dfdy(x0,y0);
    xd=@(d) x0+d*s1;
    yd=@(d) y0+d*s2;
    bd=@(d) b(xd(d),yd(d));
    d_star=fminsearch(bd,0);
    x1=xd(d_star);
    y1=yd(d_star);
    %subplot 211
    hold on
    plot(x1,y1,'o');
    pause(0.05);
    beta1=abs(x0-x1);
    beta2=abs(y0-y1);
    x0=x1;
    y0=y1;
    i=i+1;
    iter=iter+1;
    % Print only for 100 iterations
    if mod(iter,100)==0
        iter
    end
end
toc
x1
y1
b(x1,y1)
disp("Local minimum solution approach: Steepest Descent's method")
fprintf("The solution is obtained after %d iterations\n",i)
fprintf("The local minimum point : %0.5f,%0.5f\n",x1,y1)