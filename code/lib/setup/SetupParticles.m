function P = SetupParticles(S, shape, param1, param2)

    %% POSITIONS DEPENDING ON SHAPE PARAM
    if strcmp(shape,'sphere');
        R0 = param1;
        phi = unifrnd(0,2*pi, 1, S.M);
        costheta = unifrnd(-1,1, 1, S.M);
        u = unifrnd(0,1,1,S.M);
        theta = acos( costheta );
        r = R0 * u.^(1/3);

        P.x1 = .5 * S.L + r .* sin( theta ) .* cos( phi );
        P.x2 = .5 * S.L + r .* sin( theta ) .* sin( phi );
        P.x3 = .5 * S.L + r .* cos( theta );
    end
    if strcmp(shape,'uniform');
        P.x1 = unifrnd(0,S.L,1,S.M);
        P.x2 = unifrnd(0,S.L,1,S.M);
        P.x3 = unifrnd(0,S.L,1,S.M);
    end
    if strcmp(shape,'line');
        P.x1 = .5*S.L*ones(1,S.M);
        P.x2 = .5*S.L*ones(1,S.M);
        P.x3 = unifrnd(0,S.L,1,S.M);
    end

    %% VELOCITIES START AT 0
    P.u = zeros(1,S.M);
    P.v = zeros(1,S.M);
    P.w = zeros(1,S.M);

end

