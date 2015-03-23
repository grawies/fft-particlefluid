function uh = ffts(u)
uh = fftshift(fft(u));
uh(0) = 0;