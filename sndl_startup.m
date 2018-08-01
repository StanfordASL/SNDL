clear ; close all; clc;

dirs = {'function_gen_eval','kernel_features','learned_functions','learning_params','optimization',...
        'ccm_controller','data_gen','test','iLQR'};
disp('Added to path:')    
for i = 1:length(dirs)    
    addpath(genpath(strcat(pwd,filesep,dirs{i})));
    disp(dirs{i});
end
