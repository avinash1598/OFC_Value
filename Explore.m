baseDIR = '/Volumes/LaCie/recdata';

% Call the function with the base directory
listFilesInDirectory(baseDIR);

% Function to list all files in a directory and its subdirectories
function listFilesInDirectory(baseDIR)
    % Get the list of all folders and files in the base directory
    files = dir(baseDIR);
    
    % Remove '.' and '..' which represent the current and parent directories
    files = files(~ismember({files.name}, {'.', '..'}));
    
    % Loop through each file/folder in the directory
    for i = 1:length(files)
        % Get the full path of the current item
        fullPath = fullfile(baseDIR, files(i).name);
        
        % If it's a directory, recursively list the files within it
        if files(i).isdir
            disp(['Directory: ' fullPath]);
            listFilesInDirectory(fullPath); % Recursive call for subdirectory
        else
            % If it's a file, print the file name
            disp(['File: ' fullPath]);
        end
    end
end
