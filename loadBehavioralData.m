function [bhvData] = loadBehavioralData(baseDIR)
    % loadNeuronsData Load electrophysiological data based on the specified region and alignment.
    %
    %   This function loads electrophysiological data from the specified base directory.
    %   The data is filtered by the given brain region (e.g., "prefrontal cortex") and 
    %   alignment type (e.g., "stimulusOnset" or "choiceOnset"). It returns
    %   firing rates of neurons from the given region and alignment.
    %
    %   Parameters:
    %       baseDir (string or char array): The base directory where the ephys data files 
    %                                       are stored for a specific monkey for a specific session.
    %                                       Example: '/Volumes/LaCie/recdata/George00_rec29_04052021'.
    %
    %
    %   Example:
    %       baseDir = '/Volumes/LaCie/recdata/George00_rec29_04052021';
    
    % Get the list of all files and folders in the base directory
    files = dir(baseDIR);
    
    % Remove '.' and '..' which represent the current and parent directories
    files = files(~ismember({files.name}, {'.', '..'}));
    
    foundFilePath = '';

    % Loop through each file/folder and search for the correct file
    for i = 1:length(files)

        if startsWith(files(i).name, '._')
            continue;  % Skip the current iteration if file is a '._' file
        end
        
        % Get the full path of the current file/folder
        fullPath = fullfile(baseDIR, files(i).name);
        
        % If it's a directory, recursively list the files within it
        if files(i).isdir
            disp(['Directory: ', fullPath]);
            listFilesInDirectory(fullPath); % Recursive call for subdirectory
        else
            if ~isempty(foundFilePath)
                disp("Why are there multiple files in this directory! Don't know how to deal with this. Bye Bye!")
                foundFilePath = '';
                break;
            end
            % Display the found file path
            disp(['Found following file: ', fullPath]);
            foundFilePath = fullPath;
        end
    end

    % Check if a valid file was found
    if isempty(foundFilePath)
        disp('No matching file found.');
        return;
    end

    % Load the data from the found file
    bhvData = load(foundFilePath);
end
