% Checks size of an object
% from https://au.mathworks.com/matlabcentral/answers/14837-how-to-get-size-of-an-object
% not quite correct for some variable types


function s = get_size(this) 
    s = 0; 
   
    if (~isobject(this)|| isstring(this))
        % if a single variable, returns variable size
        s = whos('this').bytes;
    else
        % otherwise sum all object subcomponents
        props = properties(this);
        for ii=1:length(props) 
            currentProperty = getfield(this, char(props(ii))); 
            curr_s = get_size(currentProperty); 
            s = s + curr_s; 
        end
        % s = s * numel(this);
        % disp(s);
    end
end
