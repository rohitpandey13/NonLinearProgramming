%% Problem 1
%% f(x,y)=x^2+y^3
figure()
x=linspace(-100,100);
y=linspace(-100,100);
[xx,yy]=meshgrid(x,y);
zz=xx.^2+yy.^3;
surf(xx,yy,zz,'FaceColor','interp')
ix=find(imregionalmin(zz));
hold on
plot3(xx(ix),yy(ix),zz(ix),'r*','MarkerSize',24)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
axis square
%% g(x,y)=x^2+y^4
figure()
x=linspace(-100,100);
y=linspace(-100,100);
[xx,yy]=meshgrid(x,y);
zz=xx.^2+yy.^4;
surf(xx,yy,zz,'FaceColor','interp')
ix=find(imregionalmin(zz));
hold on
plot3(xx(ix),yy(ix),zz(ix),'r*','MarkerSize',24)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
axis square
%% problem 2
%% f(x,y)=y^2-x^2*y-2x^2*y+2x^4
figure()
x=linspace(-10,10);
y=linspace(-10,100);
[xx,yy]=meshgrid(x,y);
zz = yy.^2-(xx.^2).*yy-(2*(xx.^2)).*yy+2*(xx.^4);
surf(xx,yy,zz,'FaceColor','interp')
ix=find(imregionalmin(zz));
hold on
plot3(xx(ix),yy(ix),zz(ix),'r*','MarkerSize',24)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
axis square
shading interp
colorbar
colormap('spring')
mesh(xx,yy,zz)
%% Problem 3
%% f(x,y)=(x-2y)^2+x^4
figure()
x=linspace(-10,10);
y=linspace(-10,10);
[xx,yy]=meshgrid(x,y);
zz = (xx-2*yy).^2+(xx.^4);
surf(xx,yy,zz,'FaceColor','interp')
ix=find(imregionalmin(zz));
hold on
plot3(xx(ix),yy(ix),zz(ix),'r*','MarkerSize',24)
%% Problem 4
figure()
x=linspace(0,100);
y=linspace(-100,100);
[x1,x2]=meshgrid(x,y);
zz = 2*x1.^2+x2.^2-2*x1.*x2+2*x1.^3+x1.^4;
surf(xx,yy,zz,'FaceColor','interp')
ix=find(imregionalmin(zz));
hold on
plot3(xx(ix),yy(ix),zz(ix),'r*','MarkerSize',24)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
%% Optimal solution
fun=@(x)2*x(1).^2+x(2).^2-2*x(1).*x(2)+2*x(1).^3+x(1).^4;
x0=[0,0]
lb=[-Inf,-Inf];
ub=[Inf Inf];
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,[],[],[],[],lb,ub)
%% problem 5. a)
figure()
x=linspace(-100,100);
y=linspace(-100,100);
[x1,x2]=meshgrid(x,y);
zz = 0.5*x1.^2+x2.^2;
surf(xx,yy,zz,'FaceColor','interp')
ix=find(imregionalmin(zz));
hold on
plot3(xx(ix),yy(ix),zz(ix),'r*','MarkerSize',24)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
%% Optimal solution 
fun=@(x)0.5*x(1).^2+x(2).^2;
x0=[0,0]
A=[-2 -1;1 -1]
b=[-2; 1]
lb=[0,-Inf];
ub=[Inf Inf];
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,[],[],lb,ub)
%% Problem 5 . b)
figure()
x=linspace(0,100);
y=linspace(0,100);
[x1,x2]=meshgrid(x,y);
zz = 2*x1.^2+2*x1*x2+x2.^2-20*x1-14*x2
surf(xx,yy,zz,'FaceColor','interp')
ix=find(imregionalmin(zz));
hold on
plot3(xx(ix),yy(ix),zz(ix),'r*','MarkerSize',24)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
%% Optimal solution 
fun=@(x)2*x(1).^2+2*x(1)*x(2)+x(2).^2-20*x(1)-14*x(2);
x0=[3,4]
A=[-1 -1;2 -1]
b=[-25; 4]
lb=[0,0];
ub=[Inf Inf];
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,[],[],lb,ub)




