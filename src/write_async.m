function [] = write_async(gene_corrs, out_file)
    % This function writes gene correlation data to a file asynchronously.
    % Parameters:
    % gene_corrs - The data to be written to the file.
    % out_file - The output file path.
    
    % Try to acquire the lock file.
    lockFile = getLockFile(out_file, 120);
    % write
    writeToFile(gene_corrs, out_file);
    % Delete the lock file if it exists.
    if isfile(lockFile); delete(lockFile); end
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

function lockFile = getLockFile(out_file, timeout)
    % This function creates a lock file to prevent simultaneous writes.
    % It waits for a specified timeout period if the lock file already exists.
    %
    % Parameters:
    % out_file - The output file path.
    % timeout - The maximum time to wait for the lock file (in seconds).
    %
    % Returns:
    % lockFile - The path to the lock file.
    
    lockFile = '../results/output' + string(keyHash(out_file)) + '.lock';
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