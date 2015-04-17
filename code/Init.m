function Init(seed)

    pathStr = pwd;
    addpath(strcat(pathStr,'/lib'));
    addpath(strcat(pathStr,'/lib/fft_shifted'));
    addpath(strcat(pathStr,'/lib/setup'));
    addpath(strcat(pathStr,'/lib/deltasum'));
    addpath(strcat(pathStr,'/lib/plotting'));
    addpath(strcat(pathStr,'/scripts'));
    
    % RANDOM NUMBER GENERATOR!
    rng(seed);

end