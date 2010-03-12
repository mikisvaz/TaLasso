function [GE, ME, genes, miRNAs] = loadExpressionData(path)
  GE = dlmread([path, '/', 'geneExpression.txt']); 
  ME = dlmread([path, '/', 'mirnaExpression.txt']); 
  genes  = textread([path, '/', 'gene.txt'], '%q'); 
  miRNAs = textread([path, '/', 'mirna.txt'], '%q'); 

