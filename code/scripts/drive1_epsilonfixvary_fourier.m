Init(0);

cstring='rgbcmyk';
lstring='+osd*.';

L=4*2*pi;
epsilon = L/16;

N0 = 8;

%% LOOPING

n = 5;
m = 4000;

Nlist = zeros(n,1);
for i = n:-1:1
    % conditions
    N = N0 * 2^i;
    Nlist(n+1-i) = N; N
    % uncomment below for varying epsilon:
    %epsilon = L/N;
    
    % calculate
    [S,P,F] = CalcPreTime(L,N,epsilon,1);
    fh = fftns(F.hdeltasum);
    
    % transform to plot variables
    u=sqrt(conj(fh).*fh);
    u = u(:);
    k = sqrt(S.ksq);
    k = k(:);
    % choose random subset of m elements
    rI = randperm(numel(u));
    u = u(rI(1:m));
    k = k(rI(1:m));
    
    % sort
    [k,I] = sort(k);
    u = u(I);
    
    % plot
    linespec = strcat(cstring(i),lstring(i));
    loglog(k,u,linespec);
    hold on;
end

%% COSMETICS

ylabel('fft(deltasum)');
xlabel('|k|');
title('fourier coefficients of the delta sums');
legends = [];
for i = 1:n
    %legends = [legends; strcat('N = ',int2str(Nlist(i)))];
end
legends = strcat('N = ', int2str(Nlist));
legend(legends);
grid on;