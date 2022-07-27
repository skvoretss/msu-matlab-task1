function f = cosfunc_ob(x)
% мы хитрюги 
% и поэтому развернули функцию зеркально относительно OX
f = -cos(power(x, 2) - 4*abs(x)); 
end