%% Evolution algorithm
%% Optimal solution using fminsearch
clear all
clc
f=@(x,y) (x(1)-x(2))/((x(1)^2+5)*(x(2)^2+5))+(x(2)^2)/20000
x0=[0,0];
[XY]=fminsearch(f,x0)
%% Plotting the optimal solution
ff=@(x,y) (x-y)./((x.^2+5).*(y.^2+5))+(y.^2)/20000
fcontour(ff)
grid
axis square
hold on
plot(XY(1,1),XY(1,2),'r+','MarkerSize',10)
%% Generation 1
x_min=0;
y_min=0;
i=0;
ITMAX=1000;
    x=rand(10,1)*10-5+x_min;
    y=rand(10,1)*10-5+y_min;
while i < ITMAX
    %plot(x,y,'o');
    obj=f(x,y);
    [obj_min,I]=min(obj);
    x_min=x(I);
    y_min=y(I);
    plot(x_min,y_min,'b+')
    pause(0.05)
    x=rand(10,1)*x_min-1;
    y=rand(10,1)*y_min-1;
    i=i+1
end
x_min
y_min
ff(x_min,y_min)
