function [] = write_async(gene_corrs, out_file, ii)
    % This function writes gene correlation data to a file asynchronously.
    % Parameters:
    % gene_corrs - The data to be written to the file.
    % out_file - The output file path.
    % ii - The number of attempts to acquire the lock file (default is 1).
    
    if nargin < 3
        ii = 10;
    end
    
    if ii == 0
        writeToFile(gene_corrs, out_file);
    else
        % Try to acquire the lock file.
        lockFile = getLockFile(out_file, 120, ii);
        % Recurse
        write_async(gene_corrs, out_file, ii-1);
        % Delete the lock file if it exists.
        if isfile(lockFile); delete(lockFile); end
    end
end

function writeToFile(gene_corrs, out_file)
    % This function writes the gene correlation data to the specified file.
    % It appends the data to the file if it already exists.
    try
        writetable(gene_corrs, out_file, 'WriteMode', 'Append');
    catch ME
        % Display an error message if writing fails.
        disp('Error writing to file:');
        disp(ME.message);
    end
end

function lockFile = getLockFile(out_file, timeout, ii)
    % This function creates a lock file to prevent simultaneous writes.
    % It waits for a specified timeout period if the lock file already exists.
    %
    % Parameters:
    % out_file - The output file path.
    % timeout - The maximum time to wait for the lock file (in seconds).
    % ii - The current attempt number.
    %
    % Returns:
    % lockFile - The path to the lock file.
    
    lockFile = 'lock/output' + string(keyHash(out_file + ii)) + '.lock';
    elapsed = 0;
    while isfile(lockFile)
        % Wait for 1 second if the lock file exists.
        % pause(1);
        java.lang.Thread.sleep(1000)
        elapsed = elapsed + 1;
        if elapsed > timeout
            disp('Timeout reached, proceeding without acquiring lock.');
            break;
        end
        disp(elapsed);
    end
    
    fid = fopen(lockFile, 'w');
    if fid ~= -1
        fclose(fid);
    else
        error('Could not create lock file.');
    end
end