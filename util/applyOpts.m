%                               applyOpts.m
%
% A function to merge two structs of options into a single struct. This is
% useful for merging a struct of user-supplied options into a struct of
% default options.
%
% For a field that appears in both otherOpts and opts, if override flag is turned on, 
% this field in opts will take the value from the same field in otherOpts.
% For a field that only appears in otherOpts, it will be copied into opts.
% 
% Inputs:
%         opts(struct)       :  consists of options.
%         otherOpts(struct)  :  consists of options.
%         override(boolean)  :  If a field in opts that is also non-empty
%                               in otherOpts takes the value in otherOpts.
% Outputs:
%         opts(struct)       :  consists of eventual options
%
% 
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

function opts = applyOpts(opts, otherOpts, override)
    otherOptNames = fieldnames(otherOpts);
    
    for i = 1 : length(otherOptNames)
        optName = otherOptNames{i};
        if ~isfield(opts, optName) || override
            opts.(optName) = otherOpts.(optName);
        end
    end
end