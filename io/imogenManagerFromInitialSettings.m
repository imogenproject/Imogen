function manager = imogenManagerFromInitialSettings(initialSettings)
% This function uses an initial settings (ini) structure created by an Initializer subclass to 
% create an appropriate ImogenManager for use outside of imogen. The initial settings structure can
% either be supplied as an argument or, if no argument is specified, an attempt will be made to load
% it from the ini_settings.mat file in the current folder if such a file exists.
%
%>> initialSettings     initializer settings to use to create the manager.          struct
%<< manager             initialized imogen run manager.                             ImogenManager

    %--- Handle arguments ---%
    if nargin < 1 || isempty(initialSettings)
        if exist('ini_settings.mat','file') > 0
            S               = load('ini_settings.mat');
            fields          = fieldnames(S);
            initialSettings = S.(fields{1});
            clear('S');
        else
            error('Imogen:InputError: Unable to find initial settings variable on disk to load.');
        end
    end

    fprintf('Creating an ImogenManager from initial settings...\n');
    
    %--- Backward compatibility corrections ---%
    defaultValuesObject = Initializer(initialSettings.grid);
    fields = fieldnames(defaultValuesObject);
    for i=1:length(fields)
        if ~isfield(initialSettings,fields{i})
            initialSettings.(fields{i}) = defaultValuesObject.(fields{i});
            fprintf('Missing setting %s has been defaulted.\n',fields{i});
        end
    end
    
    manager = initialize(initialSettings);
end
