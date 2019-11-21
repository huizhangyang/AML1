function [x,reconError] = AML1(n,m,k,maxIter,A,At,xt,c,lambda,isComplex)
    function y = Afun(x,transp_flag)
       if strcmp(transp_flag,'transp')       % y = A'*x
          y = At(x);
       elseif strcmp(transp_flag,'notransp') % y = A*x
          y = A(x);
       end
    end
%c = measurement vector
%xt = true signal
%n = length of signal
%k = sparsity
%lambda = threshold
maxInnerIters=1000; %Iterations for the least squares
tol = 1e-4;


%Initial guess at x0 (as per AM paper).
J = randi([1 n],1,k);
z = zeros(n,1);
for j = J
    if(normrnd(0,1) < 0)
        z(j) = 3 + rand(1,1);
    else
        z(j) = -3 - rand(1,1);
    end
end
x=z;

for r = 1:maxIter
    %%Compute Projection of xk onto Z_c (Equation 5 in AM paper)
    magnitude = abs(A(x));
    a2 = A(x);
    for j = 1:m
 
        if magnitude(j) ~=0
            z(j) = c(j) * (a2(j)/magnitude(j));
        else
            theta(j) = normrnd(0,1);
            z(j) =  c(j)*exp(1i*theta(j));
        end
    end
     
%x(k+1) = T_{lambda}(Re(P_(Zc) (x(k)))
% Solve the least-squares problem 
%  z = \argmin ||A(xi)-z||^2.
% If A is a matrix,
%  z = inv(A)*z

% If A is a fourier transform( and measurements are not oversampled i.e. m==n),
%  z = inverse fourier transform of z  
% Use the evalc() to capture text output, thus preventing
% the conjugate gradient solver from printing to the screen.
x0 = rand(n,1);
evalc('z=lsqr(@Afun,z,tol,maxInnerIters,[],[],x0)');
z = real(z);
x = wthresh(z,'s',lambda);
%Find the best match over all trivial ambiguities

x = bestMatch(x,xt);

%Calculate reconstruction error
reconError(r) = norm(xt-x)/norm(xt);
reconError2(r) = norm(abs(A(x)) - c) / norm(c);
end
%figure(2)
%plot(1:length(xt),xt,1:length(x),x);
% plot(linspace(-10,10),linspace(-10,10), 'r-')
% hold on
% plot(x,xt, 'bo');
% ylabel('xt (True signal)')
% xlabel('xk (Reconstructed signal)')

end




