function [neuronsData] = loadNeuronsData(baseDIR, region, alignment)
    % loadNeuronsData Load electrophysiological data based on the specified region and alignment.
    %
    %   This function loads electrophysiological data from the specified base directory.
    %   The data is filtered by the given brain region (e.g., "prefrontal cortex") and 
    %   alignment type (e.g., "stimOnset" or "choiceOnset"). It returns
    %   firing rates of neurons from the given region and alignment.
    %
    %   Parameters:
    %       baseDir (string or char array): The base directory where the ephys data files 
    %                                       are stored for a specific monkey for a specific session.
    %                                       Example: '/Volumes/LaCie/recdata/George00_rec29_04052021'.
    %
    %       region (string or char array): The brain region for which the ephys data should be loaded.
    %                                      This helps filter the relevant data for that region. 
    %                                      Example: LOFC (left orbitofrontal cortex), 
    %                                               ROFC: (right orbitofrontal cortex), 
    %
    %       alignment (string or char array): The alignment type used for synchronizing the ephys data.
    %                                          This could correspond to specific events like 'stimOnset' 
    %                                          or 'choiceOnset'. It determines how the data should be 
    %                                          aligned or processed. 
    %                                          Example: 'stimOnset', 'choiceOnset'.
    %
    %   Example:
    %       baseDir = '/Volumes/LaCie/recdata/George00_rec29_04052021';
    %       region = 'LOFC';
    %       alignment = 'stimulusOnset';
    %       loadEphysData(baseDir, region, alignment);

    % Initialize alignment string based on alignment type
    alignmentString_ = '';
    
    if strcmp(alignment, 'stimOnset')
        alignmentString_ = 'pic';
    elseif strcmp(alignment, 'choiceOnset')
        alignmentString_ = 'choice';
    else
        disp('Invalid parameter for alignment. Must be one of "stimOnset" or "choiceOnset".');
        return;
    end

    % Check if the region is valid
    if ~ismember(region, {'LOFC', 'ROFC', 'LACC', 'RACC'})
        disp('Invalid parameter for region. Must be one of "LOFC", "ROFC", "LACC", "RACC".');
        return;
    end
    
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
            % If it's a file, check if it matches the region and alignment
            if contains(fullPath, region) && contains(fullPath, alignmentString_)
                % Display the found file path
                disp(['Found following file: ', fullPath]);
                foundFilePath = fullPath;
                break;
            end
        end
    end

    % Check if a valid file was found
    if isempty(foundFilePath)
        disp('No matching file found.');
        return;
    end
    
    % Load the data from the found file
    neuronsData = load(foundFilePath);
end
