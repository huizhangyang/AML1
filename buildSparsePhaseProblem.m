function [A,At,xt,c] = buildSparsePhaseProblem(m,n,k,SNR, isComplex,measType)
%n = length of signal
%k = sparsity
%isComplex = true or false
%measType = linear or square measurements
%SNR = noise

% rips first 'length' entries from a vector
rip = @(x,length) x(1:length);
A = @(x) fft([x;zeros(m-n,1)]);
At = @(x) rip(m*ifft(x),n);     % transpose of FM
J = randi([1 n],1,k);%Coordinates of non-zero elements
xt = zeros(n,1);
for j = J
    if(normrnd(0,1) < 0)
        xt(j) = 3 + rand(1,1);
    else
        xt(j) = -3 - rand(1,1);
    end
end
%Make measurement vector

if measType == "linear"
    noise = awgn(zeros(m,1),SNR);
    c = abs(A(xt)) + noise; 
else
    % Compute the phaseless noisy measurements for square measurements
    noise = awgn(zeros(m,1),SNR);
    c = abs(A(xt)).^2 + noise; 

    for j =1:length(c)
      if c(j) <0
         c(j) = 0; 
      end
    end
c = sqrt(c);
end
end

