function f=linesearch(func,dfunc,x)
% Exact line search method as described in notes 
x=x(1);
y=x(2);
fold=eval(func);
d=dfunc; % direction of descent
u=1;
delta=0.5;
alfa=1; % To initialize at 1
while u>0
    x1=x+alfa*d;
    x=x1(1);
    y=x1(2);
    fnew=eval(func);
    if fnew < fold
        %fprintf("steplength: %0.5f\n",alfa);
        u=0;
        f=alfa;
    else
        alfa=alfa*delta;
    end
end
end
