function LocalInit(seed)

    cwd = pwd;
    addpath('..');
    cd('..');
    Init(seed);
    cd(cwd);
end