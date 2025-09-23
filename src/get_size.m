% Checks size of an object
% from https://au.mathworks.com/matlabcentral/answers/14837-how-to-get-size-of-an-object
% not quite correct for some variable types


function s = get_size(obj)

    % s = whos('this').bytes;

    byteStream = getByteStreamFromArray(obj);
    s = numel(byteStream);

    disp([num2str(s) ' bytes']);
end
