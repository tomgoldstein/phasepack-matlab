%                               getExtraOpts.m
% 
% This function raises a warning for all the fields in extras struct 
% since they won't be used by the solver in estimating the unknown signal.
% It is used as a helper function in manageOptions.m and
% manageOptionsForBenchmark.m.
% 
% Inputs:
%         extras(struct)   :  consists of extra, unnecessary options 
%                             provided by the user.
% 
% 
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

function warnExtraOpts(extras)
    optNames = fieldnames(extras);
    
    for i = 1 : length(optNames)
        optName = optNames{i};
        warning('Provided option "%s" is invalid and will be ignored', optName);
    end   
end