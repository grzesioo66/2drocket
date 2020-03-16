function [ dx_dt ] = aa_rhs( x , t , u )

  dx_dt = zeros( size(x)(1) , 1 );

  deg2rad = pi/180.0;
  g = 9.81;

  dx_dt(1) = x(2);
  dx_dt(2) = -g*sin(x(1)) + u;

endfunction
