function F = simpson(y, x)
%This finction calculate integral using Simpson-method
if length(x) == length(y)
    h = (x(end) - x(1))/length(x);
    F = h/3 * (y(1) + y(end) + 2*sum(y(3:2:end-1)) + 4*sum(y(2:2:end-1)));
else
    disp('length(x) != length(y)');
end
end