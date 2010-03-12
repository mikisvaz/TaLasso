function [ Beta , Beta_sorted , sorted ] = GenMiRmodified( C , X , Z , gene , mirna )

% X = gene expression
% Z = mirna expression
% C = gene x mirna putative targets

N = size( X , 1 );
mu = mean(X,1);
alpha = 5e-4;

X = X-repmat(mu,N,1);
r = (1/N)*(sum(X.^2,1))';
Xn = X./(r*ones(1,N))';  

Beta = -alpha*(Xn*Z').*C+0.5*C;

[ I , J ] = find( Beta > 0 );
B = Beta( Beta > 0 );
[ Beta_sorted , i ] = sort( -B );
Beta_sorted = - Beta_sorted;
sorted = [ gene( I ( i ) ) , mirna( J ( i ) ) ];
