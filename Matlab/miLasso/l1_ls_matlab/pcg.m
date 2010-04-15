function [x,flag,relres,iter] = pcg_mod2(A,b,tol,maxit,M1,M2,x0,varargin)

m = size(b,1);
n = m;
    
% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
                                   % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('pcg',tol,maxit,0,flag,iter,NaN);
    end
    return
end

x = x0;

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
%r = b - A(x,varargin{:});
r = b - ((varargin{2}*((varargin{1}*x)*2))+varargin{3}.*x);
%r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    if (nargout < 2)
        itermsg('pcg',tol,maxit,0,flag,iter,relres);
    end
    return
end

normrmin = normr;                  % Norm of minimum residual
rho = 1;
stag = 0;                          % stagnation of the method
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    
    %y = M1(r,varargin{:});
    y = varargin{4}.*r;
    if ~all(isfinite(y))
        flag = 2;
        break
    end
    
    z = y;
        
    rho1 = rho;
    rho = r' * z;
    if ((rho == 0) || isinf(rho))
        flag = 4;
        break
    end
    if (ii == 1)
        p = z;
    else
        beta = rho / rho1;
        if ((beta == 0) || isinf(beta))
            flag = 4;
            break
        end
        p = z + beta * p;
    end
    
    %q = A(p,varargin{:});
    q = (varargin{2}*(varargin{1}*(p*2)))+varargin{3}.*p;
    pq = p' * q;
    if ((pq <= 0) || isinf(pq))
        flag = 4;
        break
    else
        alpha = rho / pq;
    end
    if isinf(alpha)
        flag = 4;
        break
    end
    
    % Check for stagnation of the method    
    if (norm(p)*abs(alpha) < eps*norm(x))
        stag = stag + 1;
    else
        stag = 0;
    end
    
    x = x + alpha * p;             % form new iterate
    r = r - alpha * q;
    normr = norm(r);
    normr_act = normr;
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        %r = b - A(x,varargin{:});
        r = b - ((varargin{2}*((varargin{1}*x)*2))+varargin{3}.*x);
        normr_act = norm(r);
        if (normr_act <= tolb)
            flag = 0;
            iter = ii;
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning('MATLAB:pcg:tooSmallTolerance', ...
                        strcat('Input tol may be smaller than eps*cond(A)',...
                        ' and might not be achieved by PCG\n',...
                        '         Try to use a bigger tolerance'));
                end
                flag = 3;
                iter = ii;
                break;
            end
        end
    end
    if (normr_act < normrmin)      % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    if stag >= maxstagsteps
        flag = 3;
        break;
    end
end                                % for ii = 1 : maxit

% returned solution is first with minimal residual
if (flag == 0)
    relres = normr_act / n2b;
else
    %r_comp = b - A(xmin,varargin{:});
    r_comp = b - ((varargin{2}*((varargin{1}*xmin)*2))+varargin{3}.*x);
    %r_comp = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:});
    if norm(r_comp) <= normr_act
        x = xmin;
        iter = imin;
        relres = norm(r_comp) / n2b;
    else
        iter = ii;
        relres = normr_act / n2b;
    end
end

