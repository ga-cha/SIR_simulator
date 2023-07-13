% here we describe an array for the yingqiuz atlas such that
% y_roi(index(i)) = g_roi(i)

order = [];
for i = 1:numel(g_roi)
    min = inf;
    for j = 1:numel(y_roi)
        diff = abs(g_roi(i)-y_roi(j));
        if diff < min
            min = diff;
            order(i) = j;
            fprintf('%d %d %d %d %d\n', order(i), j, g_roi(i), y_roi(j), y_roi(order(i)))
        end
    end
end
disp(order)
for i = 1:numel(g_roi)
    fprintf('%d %d %d\n', i, g_roi(i), y_roi(order(i)))
end
