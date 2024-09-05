function [] = write_async(gene_corrs, out_file)
    % TODO: Check if lock actually works. Low priority. It seems to.
    % Create file lock 
    lockFile = 'output.lock';
    lockFileID = fopen(lockFile, 'w');
    while lockFileID == -1
        pause(1);
        lockFileID = fopen(lockFile, 'w');
    end        
    
    % write to file, appending if available
    % if ~isfile(out_file)
    %     writetable(gene_corrs, out_file)
    % else
        writetable(gene_corrs, out_file, 'WriteMode', 'Append')
    % end

    % Release file lock
    fclose(lockFileID);
    if isfile(lockFile)
        delete(lockFile);
    end
end