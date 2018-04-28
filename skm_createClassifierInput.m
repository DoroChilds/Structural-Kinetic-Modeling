function skm_createClassifierInput(paramValues_train, labels_train, modelParam_names, filename, classifier, classNames, paramValues_test, labels_test)
% Create training data for the C4.5 classifier for decision trees. Optionally, input for the 
% proprietory C5.0 classifier can be created instead of C4.5. This function creates a '.names' 
% and '.data' file. Optionally, a '.test' file with separate test data can be created as well.

%% Default values:
if nargin < 6
    classNames = {'UNSTABLE', 'STABLE'};
    if nargin < 5
        classifier = 'C4.5';
    end
end    

% Concatenate parameters to one big matrix:
[X_TRAIN, FEAT_NAMES_TRAIN] = subFct_concatenateParams(paramValues_train, modelParam_names);

% Data dimensions:
[N_SAMPLES, N_FEATS]    = size(X_TRAIN);

%% Create NAMES-file
% Open file:
fclose('all');
f_temp = [filename, '.names'];
if exist(f_temp, 'file')
    delete(f_temp);
end

fid_NAMES = fopen(f_temp, 'a');

% Print class lables:
if strcmp(classifier, 'C5.0')
    fprintf(fid_NAMES, '%s.\n\n', 'Class');     % specific for C5.0 input
    fprintf(fid_NAMES, '%s:\t', 'Class');       % specific for C5.0 input
end
for i = 1:length(classNames)-1
    fprintf(fid_NAMES, '%s, ', classNames{i});
end
fprintf(fid_NAMES, '%s.\n', classNames{end});
if strcmp(classifier, 'C4.5')
    fprintf(fid_NAMES, '\n');   % specific for C4.5 input
end

% Print feature names:
for i=1:N_FEATS
    feat_temp = FEAT_NAMES_TRAIN{i};
    fprintf(fid_NAMES, '%s:\t%s\n', feat_temp, 'continuous.');
end

% Close file:
fclose(fid_NAMES);


%% Create DATA-file
% Open file:
fclose('all');
f_temp = [filename, '.data'];
if exist(f_temp, 'file')
    delete(f_temp);
end
fid_DATA = fopen(f_temp, 'a');

for i=1:N_SAMPLES
    x_temp = X_TRAIN(i,:);
    label_temp = classNames{labels_train(i)+1};
    
    if strcmp(classifier, 'C5.0')       % specific for C5.0 input
        fprintf(fid_DATA, '%s', label_temp);
        fprintf(fid_DATA, ', %2.2f', x_temp);
    elseif strcmp(classifier, 'C4.5')   % specific for C4.5 input
        fprintf(fid_DATA, '%2.2f,', x_temp);
        fprintf(fid_DATA, '%s', label_temp);
    end
    fprintf(fid_DATA, '\n');
    
end

% Close file:
fclose(fid_DATA);

% clear eig_temp f_temp fid_DATA i label_temp x_temp

%% Create TEST-file
if nargin > 6
    % Concatenate parameters to one big matrix:
    [X_TEST, FEAT_NAMES_TEST] = subFct_concatenateParams(paramValues_test, modelParam_names);
    
    % Open file:
    fclose('all');
    f_temp = [filename, '.test'];
    if exist(f_temp, 'file')
        delete(f_temp);
    end
    fid_DATA = fopen(f_temp, 'a');
    
    for i=1:size(X_TEST,1)
        x_temp = X_TEST(i,:);
        label_temp = classNames{labels_test(i)+1};
        
        if strcmp(classifier, 'C5.0')       % specific for C5.0 input
            fprintf(fid_DATA, '%s', label_temp);
            fprintf(fid_DATA, ', %2.2f', x_temp);
        elseif strcmp(classifier, 'C4.5')   % specific for C4.5 input
            fprintf(fid_DATA, '%2.2f,', x_temp);
            fprintf(fid_DATA, '%s', label_temp);
        end
        fprintf(fid_DATA, '\n');
        
    end
    
    % Close file:
    fclose(fid_DATA);
end
