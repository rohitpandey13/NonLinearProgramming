function [ x_opt, n_feval ] = nelder_mead(x,function_handle,flag )
% Idea was partially taken from Jeff Borggaard algorithms
% Used ANMS implementation of the Nelder MEAD algorithm
% Effect of convergence optimality is validated with 3<=n<100
% where n = number of variables in the given function
% Original reference to the algorithm is from the following
%% Define algorithm constants
  rho = 1;    % rho > 0
  xi  = 2;    % xi  > max(rho, 1)
  gam = 0.5;  % 0 < gam < 1
  sig = 0.5;  % 0 < sig < 1
  tolerance = 0.00000001;
  max_feval = 100;
%% Initialization
  [ temp, n_dim ] = size ( x );
  if ( temp ~= n_dim + 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'NELDER_MEAD - Fatal error!\n' );
    error('  Number of points must be = number of design variables + 1\n');
  end
  if ( nargin == 2 )
    flag = 0;
  end
  if ( flag )
    xp = linspace(-5,5,101);
    yp = xp;
    for i=1:101
      for j=1:101
        fp(j,i) = feval(function_handle,[xp(i),yp(j)]);
      end
    end
%% Plotting contour 
    figure (29)
    hold on
    contour(xp,yp,fp,linspace(0,200,25))
    if ( flag )
      plot(x(1:2,1),x(1:2,2),'r')
      plot(x(2:3,1),x(2:3,2),'r')
      plot(x([1 3],1),x([1 3],2),'r')
      %pause
      plot(x(1:2,1),x(1:2,2),'b')
      plot(x(2:3,1),x(2:3,2),'b')
      plot(x([1 3],1),x([1 3],2),'b')
    end
  end
  index = 1 : n_dim + 1;
  [f] = evaluate ( x, function_handle ); 
  n_feval = n_dim + 1;
  [ f, index ] = sort ( f );
  x = x(index,:);
%%  Nelder Mead iteration 
  converged = 0;
  diverged  = 0;
  while ( ~converged && ~diverged )
%  Compute the midpoint of the simplex opposite the worst point.
    x_bar = sum ( x(1:n_dim,:) ) / n_dim;
%  Compute the reflection point.
    x_r   = ( 1 + rho ) * x_bar ...
                - rho   * x(n_dim+1,:);
    f_r   = feval(function_handle,x_r); 
    n_feval = n_feval + 1;
%  Accept the point: 
    if ( f(1) <= f_r && f_r <= f(n_dim) )
      x(n_dim+1,:) = x_r;
      f(n_dim+1  ) = f_r;    
      if (flag)
        title('Reflection')
      end
%  Test for possible expansion.
    elseif ( f_r < f(1) )
      x_e = ( 1 + rho * xi ) * x_bar - rho * xi   * x(n_dim+1,:);
      f_e = feval(function_handle,x_e); 
      n_feval = n_feval+1;
%  Can we accept the expanded point?
      if ( f_e < f_r )
        x(n_dim+1,:) = x_e;
        f(n_dim+1  ) = f_e;
        if (flag), title('Expansion'), end
      else
        x(n_dim+1,:) = x_r;
        f(n_dim+1  ) = f_r;
        if (flag), title('Eventual reflection'), end
      end
%  Outside contraction.
    elseif ( f(n_dim) <= f_r && f_r < f(n_dim+1) )
      x_c = (1+rho*gam)*x_bar - rho*gam*x(n_dim+1,:);
      f_c = feval(function_handle,x_c); n_feval = n_feval+1;    
      if (f_c <= f_r) % accept the contracted point
        x(n_dim+1,:) = x_c;
        f(n_dim+1  ) = f_c;
        if (flag), title('Outside contraction'), end
      else
        [x,f] = shrink(x,function_handle,sig); n_feval = n_feval+n_dim;
        if (flag), title('Shrinkage'), end
      end
%  Try an inside contraction.
    else
      x_c = ( 1 - gam ) * x_bar + gam   * x(n_dim+1,:);
      f_c = feval(function_handle,x_c); 
      n_feval = n_feval+1;
%  Can we accept the contracted point?
      if (f_c < f(n_dim+1))
        x(n_dim+1,:) = x_c;
        f(n_dim+1  ) = f_c;
        if (flag), title('Inside contraction'), end
      else
        [x,f] = shrink(x,function_handle,sig); n_feval = n_feval+n_dim;
        if (flag), title('Shrinkage'), end
      end
    end  
    [ f, index ] = sort ( f );
    x = x(index,:);
%  Test for convergence
    converged = f(n_dim+1)-f(1) < tolerance;
%  Test for divergence
    diverged = ( max_feval < n_feval );   
    if (flag)
      plot(x(1:2,1),x(1:2,2),'r')
      plot(x(2:3,1),x(2:3,2),'r')
      plot(x([1 3],1),x([1 3],2),'r')
      %pause
      plot(x(1:2,1),x(1:2,2),'b')
      plot(x(2:3,1),x(2:3,2),'b')
      plot(x([1 3],1),x([1 3],2),'b')
    end
  end
  if ( 0 )
    fprintf('The best point x^* was: %d %d\n',x(1,:));
    fprintf('f(x^*) = %d\n',f(1));
  end
  x_opt = x(1,:);
  if ( diverged )
    x_opt=[Inf Inf];
    fprintf ( 1, '\n' );
    fprintf ( 1, 'NELDER_MEAD - Warning!\n' );
    fprintf ( 1, '  The maximum number of function evaluations was exceeded\n')
    fprintf ( 1, '  without convergence being achieved.\n' );
    fprintf ( "Algorithm terminated with user defined Maximum iterations\n");
  end
  return
end
%% Evaluate
function f = evaluate ( x, function_handle )
  [ temp, n_dim ] = size ( x );
  f = zeros ( 1, n_dim+1 );
  for i = 1 : n_dim + 1
    f(i) = feval(function_handle,x(i,:));
  end
  return
end
%% Shrinkage
function [ x, f ] = shrink ( x, function_handle, sig )
  [ temp, n_dim ] = size ( x );
  x1 = x(1,:);
  f(1) = feval ( function_handle, x1 );
  for i = 2 : n_dim + 1
    x(i,:) = sig * x(i,:) + ( 1.0 - sig ) * x(1,:);
    f(i) = feval ( function_handle, x(i,:) );
  end 
  return
end
