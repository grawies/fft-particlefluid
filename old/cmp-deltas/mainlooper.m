N = 8;
M = 500;
for i=1:5
    rng(13);
    N
    main(N,M,false);
    N = N*2;
end

