function [UNION, geneUnion, mirnaUnion] = loadPutativeData(paths)
  geneUnion = '';
  mirnaUnion = '';

  for i = 1:length(paths)
    path = char(paths(i));

    % Load data
    GM{i}    = sparse(dlmread([path, '/', 'putative.txt']));
    gene{i}  = textread([path, '/', 'gene.txt'] , '%q' );
    mirna{i} = textread([path, '/', 'mirna.txt'] , '%q' );

    % Find union of genes
    geneUnion  = [ geneUnion ; gene{i} ];
    mirnaUnion = [ mirnaUnion ; mirna{i} ];
  end

  geneUnion = unique( geneUnion );
  mirnaUnion = unique( mirnaUnion );

  UNION = zeros( length( geneUnion ) , length( mirnaUnion ) );

  for i = 1 : length(paths)    

    gm = sparse( zeros( length( geneUnion ) , length( mirna{i} ) ) );

    [ dummy , I ] = intersect( gene{i} , geneUnion );
    [ dummy , II ] = intersect( geneUnion , gene{i} );
    gm( II , : ) = GM{i}( I , : );

    gm2 = sparse( zeros( length( geneUnion ) , length( mirnaUnion ) ) );

    [ dummy , I ] = intersect( mirna{i} , mirnaUnion );
    [ dummy , II ] = intersect( mirnaUnion , mirna{i} );
    gm2( : , II ) = gm( : , I ); 

    GM{ i } = sparse( gm2 );  

    UNION = UNION + GM{i};
  end

  UNION( UNION > 1 ) = 1;

