clear

n = 2;  % dimension of x
m = 1;  % dimension of u

rad2deg = 180/pi;
deg2rad = pi/180;

R = eye(m,m);
Q = eye(n,n);


x = [ 0.0 ; 1.0  ];
u = [ 1.0 ];

dt = 0.01;
t  = 0.0;

for i = 1 : 300000 

  tp(i)   = t;
  yp(:,i) = x;
  up(:,i) = u;

  e = zeros(n,1);
  e(1) = x(1);
  e(2) = x(2);
  [A,B] = aa_matrices_AB( "aa_rhs" , x , t , u , n , m );
  [K,P] = lqr_m( A , B , Q , R );
  u = -K * e;

  x = aa_rk45( "aa_rhs" , x , t , dt , u );

  if  mod( i , 10 ) == 0 
    refresh;
    txt=sprintf('t=%7.3f ', t );
    plot( yp(1,:) , 'r' , yp(2,:) , 'b' );
    %plot( yp(4,:) , 100+50*up(2,:) , 'g' , yp(4,:) , 100*up(1,:) , 'g' , yp(4,:) , -gp(:) , 'r' , yp(4,:) , -yp(5,:) , 'b' , yp(4,i) , -yp(5,i) , 'bo' , 'linewidth' , 2 );
    title(txt);
  end

  t = t + dt;

end


