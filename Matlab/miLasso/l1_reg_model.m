function [ solucion , offset ] = l1_reg_model( GE , ME , GM , factor )

b = GE - median( GE , 2 )*ones( 1 , size(GE,2) );
a = ME;

solucion = sparse([],[],[],size(b,1),size( a , 1 ), 10*size(b,1));
offset = zeros(size(b,1),1);

lambda = 2*((b*a').*GM);
L = max(lambda');
lambda = L*factor;

A = a';
X = GM';
B = b';
Y = B;

A = A( : , X == 1 );

A = [ A , 10*ones( size( A , 1 ) , 1 ) ];

%calcular la norma de las muestras test
rel_tol = 0.001;
  
[ xi , status ] = l1_ls_nonneg( A , Y , lambda );

x = zeros( 1 , size( a , 1 ) );
x( 1 , X == 1 ) = xi(1:(size(A,2)-1));
%x(x<0.001)=0;
solucion = x;
offset = xi(size( A , 2 ));



