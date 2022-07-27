function C = my_add(A, B)
%This function calculates the sum of 2 matrix in classic way
if size(A) == size(B)
    [x, y] = size(A);
    C = zeros(x, y);
    for i = 1:x
        for j = 1:y
            C(i, j) = A(i, j) + B(i, j);
        end
    end
else
    disp('size(A) != size(B)');
end
end
