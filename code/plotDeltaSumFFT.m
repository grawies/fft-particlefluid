%cstring='+osd*.';

%% DATA GENERATION
specstring = {'-','--','-.',':','+'};
cstring = 'rgbcmyk';
for i = 1:5
    N = 8*2^i
    [S,f] = CalcDeltaFFT(N);
    k = sqrt(S.ksq);
    [k,I] = sort(k(:));
    f = f(:);
    f = f(I);
    linespec = strcat('k',specstring{i});
    plot(k(:),abs(f(:))/(N^3),cstring(i));
    hold on;
    drawnow;
end

%% COSMETICS

ylabel('fft(deltasum)');
xlabel('|k|');
title('fourier coefficients of the delta sums');
legend('one','two','three','four');