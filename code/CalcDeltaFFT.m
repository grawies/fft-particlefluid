function [S,transf] = CalcDeltaFFT(N)

    Init(0);

    L=4*2*pi;
    epsilon = L/16;
    S = SetupWorld(L, 0.4,... 
                    N, 1000,... % N,M
                    0.01, 1000,... % dt, nsteps
                    epsilon);

    P = SetupParticles(S, 'sphere', 0.5, 0.0);

    F = SetupForces(S);

    F.hdeltasum = CalcDeltaSum(S, P, F);

    transf = fftns(F.hdeltasum);
end