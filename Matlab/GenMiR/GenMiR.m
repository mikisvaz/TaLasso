% /filename GenMiR.m 
% /description Implementation of the GenMiR++ model & Variational Bayes algorithm for learning miRNA-target interactions from sequence and expression data
% /vars   X           NxT matrix of mRNA expression data
%         Z           MxT matrix of miRNA expression data
%         C           MxN matrix of putative regulatory relations bewteen miRNAs and mRNAs
%         Beta        MxN matrix of variational target selection parameters
%         Nu          Mx1 variational parameter vector for regulatory weights distribution
%		  Omega,Phi   TxT mean and variance diagonal parameter matrices for variational tissue scaling coefficients distribution
%		  a,s         Hyperparameters for regulatory weights and tissue-scaling coefficients
%         mu,Sigma    Tx1 background transcription rate and TxT diagonal covariance matrix
%         Pi          Base target selection probability  
%         glbl, mlbl  mRNA and miRNA labels, respectively
%
% /params infile   Name of MATLAB file containing mRNA matrix X, 
%                  miRNA matrix Z, a paired list of predicted miRNA-target
%                  interactions in form of binary adjacency matrix C,
%                  corresponding mRNA and miRNA labels in glbl, mlbl
%         display  Boolean flag for displaying data
%         num_iter Maximum number of iterations (default 100)
%         tol      Stopping criterion  (default 1e-5)
%         filename Name of file to which output is written (w/out extension)
%
% /version 3.0, January 2007
% /author Jim Huang, PSI Group, University of Toronto
%
% /In order to be used for this web-page, data loading and displays have been changed, by
% /Ander Muniategui, Dpto. de Electronica y Comunicaciones, CEIT (Centro de
% /Estudios e Investigaciones Tecnicas de Guipuzcoa), Universidad de Navarra
% $Id$

function [Beta Score] = GenMiR(X,Z,C)

num_iter = 100;
tol = 1e-8;

%Load data;
N = size( X , 1 );
M = size( Z , 1 );
T = size( X , 2 );


[row col] = find(C);
K = length(row);

%Initialize parameters
I = eye(T);
s = 5e-3;
mu = mean(X,1);
alpha = 5e-4;
Nu = alpha*ones(M,1);
Pi = 0.5;
Beta = Pi*C;
Omega = I;
Phi = s*I;
Beta_mat = [];

%Compute D(q||p)
Y0 = (Beta*(Z.*repmat(Nu,1,T)));
Y = Y0*Omega;
R = X-repmat(mu,N,1);
V0 = zeros(T,T);

for kk=1:K

	ii = row(kk);
	jj = col(kk);
	zj = Z(jj,:)';
	lj = Nu(jj);
	bij = Beta(ii,jj);
	V0 = V0 + (lj^2)*(2*bij-bij^2)*diag(zj.^2);	
    
end;

U = (1/N)*diag(sum(R.^2 + 2*R.*Y,1));
V = (1/N)*(Omega.^2+Phi)*V0;
W = (1/N)*(Omega.^2+Phi)*diag(sum(Y0.^2,1));
Sigma = diag(diag(U+V+W));

[Eval E dE] = GenMiR_evalE(1e30,[],C, Beta, Pi, Z, X, mu, Nu, alpha, Sigma, Omega, Phi, s, U, V, W, M, N, T, I);
t = 0;
n_iters = 1;

%Start GenMiR++

tic
while (dE > tol & t < num_iter)
		
	%Variational Bayes E-step: Infer miRNA-target interactions	
	t=t+1
 	Beta = GenMiR_VBEStep(C, Beta, Pi, Z, X, mu, Nu, Sigma, Omega, Phi, M, N, T);

	%Variational Bayes M-step
 	[Pi Nu alpha Omega Phi mu Sigma U V W] = GenMiR_VBMStep(C, Beta, Pi, Z, X, mu, Nu, alpha, Sigma, Omega, Phi, s, M, N, T, I);
			
	%Evaluate bound on likelihood
	[Eval E dE] = GenMiR_evalE(E,Eval,C, Beta, Pi, Z, X, mu, Nu, alpha, Sigma, Omega, Phi, s, U, V, W, M, N, T, I);

end;


%Score miRNA-target interactions
Score = zeros(nnz(C),1);
for kk=1:K

	ii = row(kk);
	jj = col(kk);
	
	bij = Beta(ii,jj);
	Score(kk) = log10((bij+eps)/(1-bij+eps));
    
end;
