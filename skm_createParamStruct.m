function paramIntervals = skm_createParamStruct(N, kinetics)
% skm_paramIntervals = skm_createParamStruct(N, kinetics)
% Assistant function creating a template struct that defines network positions and sampling 
% intervals of the model parameters. This struct then serves as an input argument for the 
% skm-function. Parameters can adapted by the user by manually modifying the corresponding 
% fields in the struct.

if nargin<2
    kinetics = 'enzymatic_irrev';
end
paramIntervals = subFct_checkParamIntervals(struct, N, kinetics);