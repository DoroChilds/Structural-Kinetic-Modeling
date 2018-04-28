function options = skm_createOptionsStruct()
% skm_options = skm_createOptionsStruct()
% Assistant function creating a template struct that contains all default settings for the 
% toolbox. This struct then serves as an input argument for the skm-function. Settings can 
% be adapted by the user by manually modifying the corresponding fields in the struct.


% Assign default options:
options = subFct_checkOptions(struct);
