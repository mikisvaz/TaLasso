function [ gene_X_mirna , GE , ME , DDBB_g , DDBB_m ] = GMdataIntersect( gm , DDBB_g , DDBB_m , GE_g , ME_m , GE , ME )

  [ DDBB_g , I1 , I2 ] = intersect( DDBB_g , GE_g );   
  
  GE = GE( I2 , : );   
  gene_X_mirna = gm( I1 , : );      
  
  [ Mirnas , I1 , I2 ] = intersect( ME_m , DDBB_m );
  
  DDBB_m = Mirnas;
  gene_X_mirna = gene_X_mirna( : , I2 );
  ME = ME( I1 , : );
  
  %eliminar los mirnas normalizadores con expresion = 0
  ME_m = ME_m( sum( ME' == 0 ) < size( ME , 2 ) );
  ME = ME( sum( ME' == 0 ) < size( ME , 2 ) , : );
  gene_X_mirna = gene_X_mirna( : , sum( ME' == 0 ) < size( ME , 2 ) );
  DDBB_m = DDBB_m(sum( ME' == 0 ) < size( ME , 2 ));
  
  %eliminar los genes que no tienen targets
  GE = GE( sum( gene_X_mirna' ) ~= 0 , : );
  XX = GE - mean( GE , 2 )*ones( 1 , size( GE , 2 ) );
  DDBB_g = DDBB_g( sum( gene_X_mirna' ) ~= 0 , 1 );
  gene_X_mirna = gene_X_mirna( sum( gene_X_mirna' ) ~= 0 , : ); 
  
