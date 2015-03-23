function uh = fftns(u)
uh = fftshift(fftn(u));
uh(:,1,1) = zeros(length(uh(:,1,1)),1,1);
uh(1,:,1) = zeros(1,length(uh(1,:,1)),1);
uh(1,1,:) = zeros(1,1,length(uh(1,1,:)));

