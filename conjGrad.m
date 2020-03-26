function [x,numIter]=conjGrad(func,x,b,epsilon)
if (nargin < 3 || nargin > 4)
    disp('Error');
    return;
end
if nargin==3 
     epsilon=1.0e-8;
end
n=length(b);
r=b-feval(func,x);
s=r;
for numIter=1:n
    u=feval(func,s);
    alpha=dot(s,r)/(u.*s);
    x=x+alpha*s;
    r=b-feval(func,x);
    if sqrt(dot(r,r))<epsilon
        return
    else
        beta=-dot(r,u)/(u.*s);
        s=r+beta*s;
    end
end
error('Too many iterations')
    
