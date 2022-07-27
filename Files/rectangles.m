function f = rectangles(y, x)
%This finction calculate integral using rectangles-method
if length(x) == length(y)
    h = (x(end) - x(1))/length(x);
    f = h * sum(y(1:end-1));
else
    disp('length(x) != length(y)');
end
end