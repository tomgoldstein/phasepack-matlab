%                               getExtraOpts.m
% 
% This function creates and outputs a struct that consists of all the
% options in opts but not in otherOpts.
% It is used as a helper function in manageOptions.m and
% manageOptionsForBenchmark.m
%
% Inputs:
%         opts(struct)       :  consists of options.
%         otherOpts(struct)  :  consists of options.
% Outputs:
%         extras(struct)     :  consists of extral options appearing in
%                               opts but not in otherOpts.
%
% 
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

function extras = getExtraOpts(opts, otherOpts)
    extras = struct;
    optNames = fieldnames(opts);
    
    for i = 1 : length(optNames)
        optName = optNames{i};
        if ~isfield(otherOpts, optName)
            extras.(optName) = opts.(optName);
        end
    end
end