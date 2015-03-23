function uh = fft2s(u)
uh = fftshift(fft2(u));
uh(:,1) = zeros(length(uh(:,1)),1);
uh(1,:) = zeros(1,length(uh(:,1)));

