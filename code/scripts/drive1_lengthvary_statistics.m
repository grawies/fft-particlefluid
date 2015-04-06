Init(0);

cstring='rgbcmyk';
lstring='+osd*.';

L0 = pi;
epsilon = L/16;

N0 = 8;

%% LOOPING

n = 5;

Llist = zeros(n,1);
for i = n:-1:1
    % conditions
    N = N0 * 2^i;
    L = L0 * 2^i;
    Llist(n+1-i) = L; N
    
    % calculate
    [S,P,F] = CalcPreTime(L,N,epsilon,3);
    
    tantheta = sqrt(P.u.^2 + P.v.^2) ./ P.w;
    theta = atan(tantheta);
    
    
    cdfplot(theta);
    hold on;
end

%% COSMETICS

ylabel('cdf');
xlabel('angle to z-axis');
title('Distribution of angles to Z-axis');
for i = 1:n
end
legends = strcat('L = ', num2str(Llist));
legend(legends);
grid on;