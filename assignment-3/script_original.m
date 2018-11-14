% CGtest8.m is a script for comparing {cg, symmlq, minres}
% on sparse matrix Boeing/bcsstm34, n=588, nnz=24270).
% See http://faculty.cse.tamu.edu/davis/welcome.html (Tim Davis).
% The matrix A is from a structural problem.
% It is symmetric indefinite with lambda_min = -2.6830.
%
% 09 Apr 2017: Problem Boeing/bcsstm34 used for Homework 3.
% 22 Apr 2018: Removed pcg and minres on C^2x = Cb.
%----------------------------------------------------------

load bcsstm34.mat;   % lambda(min) = -2.6830, lambda(max) = 6.8331
A     = Problem.A;   % Save original matrix A
condest(A)
[n,n] = size(A);
x     = 1./(1:n)';

%----------------------------------------------------------
sigma1 = 2.7;
B      = A + sigma1*speye(n);   condB = condest(B);
b      = B*x;
tol    = 1e-9;   % Not highly accurate
maxit  = 1000;

[xC,flagC,relresC,iterC,resvecC]       = pcg   (B,b,tol,maxit);
[xL,flagL,relresL,iterL,resvecL]       = symmlq(B,b,tol,maxit);
[xM,flagM,relresM,iterM,resvecM]       = minres(B,b,tol,maxit);
[xS,flagS,relresS,iterS,resvecS,lsvec] = lsqr  (B,b,tol,maxit);

errC  = norm(xC-x,inf);   % The inf-norm is best for large vectors
errL  = norm(xL-x,inf);
errM  = norm(xM-x,inf);
errS  = norm(xS-x,inf);

fprintf('\nPOS-DEFINITE B = A + sigma1*I,')
fprintf('  sigma1 =%5.2f,  condest(B) = %8.1e\n\n', sigma1, condB)
fprintf('                flag  iter   relres    error\n')
fprintf(' CG       Bx = b%4g %5g %8.1e %8.1e  b\n', flagC,iterC,relresC,errC)
fprintf(' SYMMLQ   Bx = b%4g %5g %8.1e %8.1e  r\n', flagL,iterL,relresL,errL)
fprintf(' MINRES   Bx = b%4g %5g %8.1e %8.1e  g\n', flagM,iterM,relresM,errM)
fprintf(' LSQR     Bx = b%4g %5g %8.1e %8.1e  m\n', flagS,iterS,relresS,errS)

figure(1)
hold off;  plot(log10(resvecL),'r-')
hold  on;  plot(log10(resvecC),'b-')
hold  on;  plot(log10(resvecM),'g-')
hold  on;  plot(log10(lsvec)  ,'m-')
xlabel('Number of iterations') % x-axis label
ylabel('Log base 10 of norm of residuals') % y-axis label
legend('symmlq','pcg', 'minres', 'lsqr')
title('Figure 1: Log base 10 of norm of residuals vs number of iterations')

%----------------------------------------------------------
sigma2 = 0.5;
C      = A + sigma2*speye(n);   condC  = condest(C);
b      = C*x;    Cfun  = @(x) C*x;       % Treat  C  as a function
%b2    = C*b;    Cfun2 = @(x) C*(C*x);   % Treat C*C as a function

[xC,flagC,relresC,iterC,resvecC]       = pcg   (Cfun ,b ,tol,maxit);
[xL,flagL,relresL,iterL,resvecL]       = symmlq(Cfun ,b, tol,maxit);
[xM,flagM,relresM,iterM,resvecM]       = minres(Cfun ,b ,tol,maxit);
[xS,flagS,relresS,iterS,resvecS,lsvec] = lsqr  (C    ,b ,tol,maxit);

errC  = norm(xC-x,inf);
errL  = norm(xL-x,inf);
errM  = norm(xM-x,inf);
errS  = norm(xS-x,inf);

fprintf('\n  INDEFINITE C = A + sigma2*I,')
fprintf('  sigma2 =%5.2f,  condest(C) = %8.1e\n\n', sigma2, condC)
fprintf('                flag  iter   relres    error\n')
fprintf(' CG       Cx = b%4g %5g %8.1e %8.1e  b\n', flagC,iterC,relresC,errC)
fprintf(' SYMMLQ   Cx = b%4g %5g %8.1e %8.1e  r\n', flagL,iterL,relresL,errL)
fprintf(' MINRES   Cx = b%4g %5g %8.1e %8.1e  g\n', flagM,iterM,relresM,errM)
fprintf(' LSQR     Cx = b%4g %5g %8.1e %8.1e  m\n', flagS,iterS,relresS,errS)

figure(2)
hold off;  plot(log10(resvecC),'b-')
hold  on;  plot(log10(resvecL),'r-')
hold  on;  plot(log10(resvecM),'g-')
hold  on;  plot(log10(lsvec)  ,'m-')

xlabel('Number of iterations') % x-axis label
ylabel('Log base 10 of norm of residuals') % y-axis label
legend('pcg', 'symmlq', 'minres', 'lsqr')
title('Figure 3: Log base 10 of norm of residuals vs number of iterations')

%----------------------------------------------------------
% Plot the eigenvalues of C.
%----------------------------------------------------------
lambda = eig(full(C));
figure(3)
hold off;    plot(lambda,'b.')
xlabel('Eigenvalue number');   ylabel('\lambda(C)');
title('Eigenvalues of C');

% Show if the eigenvalues are clustered.
figure(4)
hold off;    plot(lambda,250*ones(n,1),'b.')
hold on

y1    = -3;      yn    =  7;
step  =  0.25;   nbar  = (yn - y1)/step + 1;
y     = zeros(nbar,1);
nlam  = zeros(nbar,1);

for i = 1:nbar
  y2      = y1 + step;
  nlam(i) = length( find(lambda>y1 & lambda<=y2) );
  y(i)    = y1 + 0.5*step;
  y1      = y2;
end

bar( y, nlam )
xlabel('\lambda(C)');   ylabel('No. of \lambda(C)');
title('Distribution of eigenvalues of C');

[V, D] = eig(full(C));
c = V\b;
figure(5)
plot(sort(c, 1, 'descend'))
xlabel('Index for length of c') % x-axis label
ylabel('Coefficients of c') % y-axis label
title('Figure 7: Coefficients of c such that Vc = b')


