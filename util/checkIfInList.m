%                               applyOpts.m
%
% This function checks if an element is in a given list. If it is not, 
% an error will be raised. It is used in any place where a value is 
% expected to be chosen from a list of options
% (e.g. It is used to check opts.indexChoice in the function 
% validateOptions in solveCoordinateDescent.m).
% 
% Inputs:
%         fieldName(string)       :  The name of the field to be checked.
%         element(string)         :  The name of the element to be checked.
%         list(cell array)        :  The list that consists of all the 
%                                    valid names.
%         customizedMsg(string)   :  A customized message put at the
%                                    beginning of the error message.
%  
% 
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

function checkIfInList(fieldName,element,list,customizedMsg)
    if ~exist('customizedMsg','var')
        customizedMsg = '';
    end
    if ~any(strcmp(list,lower(element)))
        % Use num2str in case element is given as a number
        error('%s %s is invalid for field %s\n',customizedMsg,num2str(element),fieldName);
    end
end