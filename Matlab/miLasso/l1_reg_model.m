function [ solucion ] = l1_reg_model( GE , ME , GM , factor )

  b = GE - median( GE , 2 )*ones( 1 , size(GE,2) );
  a = ME;
  
  solucion = zeros( 1 , size( a , 1 ) );
  offset = 0;
  
  lambda = 2*((b*a').*GM);
  L = max(lambda');
  lambda = L*factor;
  
  for gene = 1:size(b,1)

      gene

      A = a';
      X = GM';
      B = b';
      Y = B;

      %caso general
      X = X( : , gene );
      A = A( : , X == 1 );
      B = B( : , gene );
      Y = Y( : , gene );

      A = [ A , 10*ones( size( A , 1 ) , 1 ) ];

      %calcular la norma de las muestras test
      rel_tol = 0.001;

      [ xi , status ] = l1_ls_nonneg( A , Y , lambda(gene) , rel_tol );     

      x = zeros( 1 , size( a , 1 ) );
      x( 1 , X == 1 ) = xi(1:(size(A,2)-1));
      solucion = [ solucion ; x ];
      offset = [ offset , xi(size( A , 2 )) ];

  end

  solucion = solucion( 2 : end , : );
  offset = offset( 2 : end );

