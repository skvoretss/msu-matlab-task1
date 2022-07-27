%#ok<*NOPTS>
%% Script 1
clc;

a = 0.1;
b = 6.1;
n = 100;
%n = 0.001;
%x = a:n:b;
x = linspace(a, b, n);
y = cos(power(x, 2) - 4*abs(x));

[x1, y1] = fminbnd(@cosfunc, a, b)  
[c, d] = fminbnd(@cosfunc_ob, a, b);
d = abs(d); 
x2 = c 
y2 = d
plot(x, y);
text(x1, y1, 'min'); 
text(x2, y2, 'max'); 
grid on;

%% Script 2
clc;

%2.0:
prompt = 'Input n:';
n = input(prompt);
if isnan(n)
    disp('���� � n - ��� NaN')
elseif n == Inf
    disp('���� � n - ��� Inf')
else
    if n == 2
        disp('��� ����� �������')
    elseif rem(power(n + 1, n - 1) - 1, n) == 0
        disp('��� ����� �������')
    else
        disp('��� ����� ���������')
    end

    %2.1:
    7:14:n

    %2.2:
    a = 2:n+1;
    a = a';
    d = ones(1,n+1);
    % ������-������� nx1 �������� �� ������-������ 1xn
    a * d

    %2.3: 
    b = 1:(n+1)^2;
    % �������������, ������ ��� reshape �������� ����������� 
    B = transpose(reshape(b, [n+1, n+1])) 
    % ������������� B, ����� ��������� �������� reshape
    C = reshape(transpose(B), [1, (n+1)^2])
    D = B(:, n:n+1)
end

%% Script 3
clc;
A = round(-10.5 + (10.5+10.5)*rand(5, 8))
max(diag(A))

%���������������� ����, ������ ��� �������� ����� ������� �� ��������, �
%��� ����� ������
s = sum((A'))
p = prod((A'))
if length(find(s)) ~= length(s)
    disp('������� � �������� �� ����� ���� �������, ���� ������� �� ����');
else
    vec = p./s
    d_max = max(vec)
    d_min = min(vec)
end
N = sortrows(A, 'descend')

%% Script 4 
clc;
prompt = 'Input n:';
n = input(prompt);
prompt = 'Input m:';
m = input(prompt);

if isnan(n)
    disp('���� � n - ��� NaN')
elseif n == Inf
    disp('���� � n - ��� Inf')
elseif isnan(m)
    disp('���� � m - ��� NaN')
elseif m == Inf
    disp('���� � m - ��� Inf')
else
    if rem(m, 2) == 1
        disp('������ ��������� ��������� ������� G')
    else
        %����� ������ �������� �� 1 �� n*m
        b = 1:n*m; 

        %������ ������ � ���, ����� ����� ���������� �� ��������
        %m � n ���������� �������, ����� ��������� �������� ������� A
        %(����� ����������� � ����� ���� �� ��������, � �� ��������)
        %���������������, ����� ��� ������� �� ���� �����
        A = transpose(reshape(b, [m, n]))

        %����� ��� ������ � ���������� ����� 1 (�.�. ����������� � g1 �����
        %��������� ������� �)
        g1 = A(1:n, 1:2:m);

        %����� ��� ������ �������� �� ������ �������
        g2 = A(2:2:n, 2:2:m);

        %"������������ ��������" ������� � �������� �� ������ ��� ��-�� g2
        g1(2:2:end, :) = g2

        B = A(2:2:n, 1:2:m)
        R = A(1:2:n, 2:2:m)
    end
end
%% Script 5
clc;
prompt = 'Input n here:';
n = input(prompt);
prompt = 'Input m here:';
m = input(prompt);
b1 = 1:n % �������� ������� �
b2 = 1:m % �������� ������� y
% �������� �������� ������� y n ��� � ����������� n
Z = repmat(b2, n);

%����� ������ ������ ������, �.�. � ��� ��� � ����������� n 
% � ������������� � �������
Y = transpose(Z(1, :));

%�������� �������� ������� x m ��� � ����������� m
%������������� b1, ����� ������� ��� �� ��������
Z = repmat(transpose(b1), m);

%����� ������ n �����(�.�. ����������� �� m)
%��������� � ������-������� - ������������� ��-�� reshape(����� ������
%�������� � �������� ���������, �.�. ������ ��� � �������)
X = reshape(transpose(Z(1:n, :)), [m*n, 1]);
A = [X Y]

%% Script 6 
clc;
%��� ����� ������ �������, ��� �����
prompt = 'Input n:';
n = input(prompt);
if isnan(n)
    disp('���� � n - ��� NaN')
elseif n == Inf
    disp('���� � n - ��� Inf')
else
    
    b = reshape(1:3*n, [3, n])
    %b = zeros(3, n)+.1;
    x = b(1, :);
    y = b(2, :);
    z = b(3, :);

    %Aij = x(yizj - ziyj) + y(zixj - xizj) + z(xiyj - yixj)

    x1 = z'.*y; % ziyj
    x2 = x1'; % yizj
    X = x2 - x1;


    y1 = x'.*z; % xizj
    y2 = y1'; % zixj
    Y = y2 - y1;


    z1 = y'.*z; % yixj
    z2 = z1'; % xiyj
    Z = z2 - z1;

    A = sqrt(X.*X + Y.*Y + Z.*Z)
    
end

%% Script 7.1
clc;
prompt = 'Input n:';
n = input(prompt);
prompt = 'Input m:';
m = input(prompt);
if isnan(n)
    disp('���� � n - ��� NaN')
elseif n == Inf
    disp('���� � n - ��� Inf')
elseif isnan(m)
    disp('���� � m - ��� NaN')
elseif m == Inf
    disp('���� � m - ��� Inf')
else
    a = -5 + (5+5)*rand(1, n)
    b = 1:2:m*2
    a_min = min(a);
    a_max = max(a);
    b_min = min(b);
    b_max = max(b);
    c1 = a_max - b_min;
    c2 = b_max - a_min;
    if c1*c2 < 0
        if c1 + c2 < 0
            res = min(c1, c2) *(-1);
        else
            res = max(c1, c2);
        end
    elseif (c1 < 0) && (c2 < 0) 
        if c1 - c2 < 0
            res = c1 *(-1);
        else
            res = c2 *(-1);
        end

    elseif (c1 > 0) && (c2 > 0) 
        if c1 - c2 < 0
            res = c2;
        else
            res = c1;
        end
    elseif c1 == 0
        if c2 < 0
            res = c2 * (-1);
        else
            res = c2;
        end
    elseif c2 == 0
        if c1 < 0
            res = c1 * (-1);
        else
            res = c1;
        end 
    end
    disp(res)
end


%% Script 7.2
clc;
prompt = 'Input n:';
n = input(prompt);
prompt = 'Input m:';
m = input(prompt);
if isnan(n)
    disp('���� � n - ��� NaN')
elseif n == Inf
    disp('���� � n - ��� Inf')
elseif isnan(m)
    disp('���� � m - ��� NaN')
elseif m == Inf
    disp('���� � m - ��� Inf')
else
    a = -5 + (5+5)*rand(1, n)
    b = 1:2:m*2
    mx = max([max(b) - min(a), max(a) - min(b)]) 
end
%% Script 8
clc;
prompt = 'Input n:';
n = input(prompt); % ���������� ��������
prompt = 'Input k:'; % ���������� ���������
k = input(prompt);

if isnan(n)
    disp('���� � n - ��� NaN')
elseif n == Inf
    disp('���� � n - ��� Inf')
elseif isnan(k)
    disp('���� � k - ��� NaN')
elseif k == Inf
    disp('���� � k - ��� Inf')
else
    %   < - - - - k - - - - >
    %   ?
    %   |
    % n |
    %   |
    %   ?
    b = reshape(1:k*n, [n, k])
    % �������������, ����� ��� ����� �������� repmat
    c = b';
    % ��������� � ������e � ������������� �������
    c = c(:);
    c = c';
    %x1 x2 x3 x4 x5...y1 y2 y3 y4 y5... - ���������� 1-�� � 2-�� �������
    %x1 x2 x3 x4 x5...y1 y2 y3 y4 y5... - ���������� 1-�� � 2-�� �������
    d1 = repmat(c, n);
    %����� ������ n �����, ������ ��� � ����� ��� ������
    d1 = d1(1:n, 1:n*k);
    d2 = repmat(b, n);
    %������ ��� � �� �����
    %x1 x2 x3 x4 x5... x1 x2 x3 x4 x5... - ���������� 1-�� �������
    %y1 y2 y3 y4 y5... y1 y2 y3 y4 y5... - ���������� 2-�� �������
    d2 = d2(1:n, 1:n*k);
    %�������� ���������
    d = d1 - d2;
    %�������� �������� � �������
    d = d.^2;
    %������������� ��� reshape
    d = d';
    %������ reshape, ����� ����� ������ �����
    d = reshape(d, [k, n*n]);
    d = sum(d);
    d = reshape(d, [n, n]);
    disp(sqrt(d))
end
%% Script 9
clc;
prompt = 'Input n:';
n = input(prompt);

if isnan(n)
    disp('���� � n - ��� NaN')
elseif n == Inf
    disp('���� � n - ��� Inf')
else 
    A = randi(5, n)
    B = randi(10, n)
    
    d = 20;
    time = zeros(d, 1);
    disp('���� ��������:');
    tic;
    for w = 1:d
        C = my_add(A, B);
        time(w) = toc;
    end
    disp(median(time));
    
    disp('�������� ���� ������� �������:');
    tic
    for w = 1:d
        C = A + B;
        time(w) = toc;
    end
    disp(median(time))
    
    z = 1;
    T1 = zeros(1,n*2);
    for q = 1:n*2
        A = 1:100*z;
        B = 1:100*z;
        tic
        for i = 1:d         
            C = A + B;
            time(i) = toc;
        end
        T1(q) = median(time);
        z = z*2;
    end

    z = 1;
    T2 = zeros(1,n*2);
    for q = 1:n*2
        A = 1:100*z;
        B = 1:100*z;
        tic
        for i = 1:d         
            C = my_add(A, B);
            time(i) = toc;
        end
        T2(q)= median(time);
        z = z*2;
    end

    plot(1:n*2, T2, 'r'); % ���� ��������
    hold on; 
    plot(1:n*2, T1, 'g'); % �������� ���� ������� �������
    xlabel('���');
    ylabel('������ �');
    legend('���� ��������', '����������� ��������');
    hold off;
    %grid on;
end

%% Script 10 
clc;
a = [1, 2, 3, 4]
%a = [1, 2, 3, 4, 5]
%a = [1, 2, 3, 3, 2, 1]
%a = [1, 2, 3, 4, 3, 2, 1]
%a = [7]
a1 = flip(a)
if isequal(a, a1)
    disp('Vector is symmetrical')
else
    disp('Vector is not symmetrical')
end

%% Script 11
clc;
prompt = 'Input n:';
n = input(prompt);
prompt = 'Input a:';
a = input(prompt);
prompt = 'Input b:';
b = input(prompt);
if isnan(n)
    disp('���� � n - ��� NaN')
elseif n == Inf
    disp('���� � n - ��� Inf')
elseif isnan(a)
    disp('���� � a - ��� NaN')
elseif a == Inf
    disp('���� � a - ��� Inf') 
elseif isnan(b)
    disp('���� � b - ��� NaN')
elseif b == Inf
    disp('���� � b - ��� Inf') 
else
    % �������� �� 0, � ����� ��������� ����� � -1, ����� ������� ������ � 
    % ������� [0, a], � �� � �������� (0, a)
    x = realmin*(-1) + a.*rand(1, n)
    c = a/2*b
    % � test ��������� ��� �������, � ������� ��-�� �� x ������ b
    test = x > b
    disp('res := percentage of elements > b');
    res = sum(test)/n
    if res >= c
        disp('res >= a/2b');
    else
        disp('res < a/2b');
    end
end

%% Script 12
clc;
v_rec = zeros(1, 30);
v_sim = zeros(1, 30);
v_trap = zeros(1, 30);


d = 30;
time = zeros(30, 1);
disp('Rectangles time:');
tic
for q = 1:d
    j = 1;
    for i = 0:1/10:3-1/10
        x = i:1/100:i+1/10;
        y = exp_square(x);
        v_rec(j) = rectangles(y, x);
        j = j+1;
    end
    time(q) = toc;
end
disp(median(time));


disp('Simson time:');
tic
for q = 1:d
    j = 1;
    for i = 0:1/10:3-1/10
        x = i:1/100:i+1/10;
        y = exp_square(x);
        v_sim(j) = simpson(y, x);
        j = j+1;
    end
    time(q) = toc;
end
disp(median(time));


disp('Trapz time:');
tic
for q = 1:d
    j = 1;
    for i = 0:1/10:3-1/10
        x = i:1/100:i+1/10;
        y = exp_square(x);
        v_trap(j) = trapz(x, y);
        j = j+1;
    end
    time(q) = toc;
end
disp(median(time));

x = 0:1/10:3-1/10;
plot(x, v_rec, 'g'); 
hold on;
plot(x, v_sim, 'r');
xlabel('���');
ylabel('��������');
legend('����� ���������������', '����� ��������');
hold off;

j = 1;
for i = 0:1/10:3-1/10
    
    x_h = i:1/50:i+1/10;
    x_h_2 = i:1/100:i+1/10;
    
    y_h = exp_square(x_h);
    y_h_2 = exp_square(x_h_2);
    
    v_sim(j) = simpson(y_h, x_h) - simpson(y_h_2, x_h_2);
    v_rec(j) = rectangles(y_h, x_h) - rectangles(y_h_2, x_h_2);
    v_trap(j) = trapz(x_h, y_h) - trapz(x_h_2, y_h_2);
    j = j+1;
end

figure
x = 0:1/10:3-1/10;
plot(x, v_rec, 'g'); 
hold on;
plot(x, v_sim, 'r');
plot(x, v_trap, 'b');
xlabel('���');
ylabel('�����������');
legend('����� ���������������', '����� ��������', '����� ��������');
hold off;

%% Script 13
clc;
y = @(x) cos(x);
y_derivative = @(x) sin(x)*(-1);
y_derivative_central = @(x, h) (y(x + h) - y(x - h)) / (2*h);
y_derivative_right = @(x, h) (y(x + h) - y(x)) / h;

x = pi / 4;
%h = logspace(-7, 0); % ������� ������ ���� � ����� 
%h = logspace(-50, 0); % ������� ����
h = logspace(-5, 0, 50);

y_diff_right = zeros(1, 50);
y_diff_central = zeros(1, 50);

for i = 1:50
    y_der = y_derivative(x);
    y_diff_right(i) = abs(y_derivative_right(x, h(i)) - y_der);
    y_diff_central(i) = abs(y_derivative_central(x, h(i)) - y_der);
end

loglog(h, y_diff_right, 'b'); 
hold on;
loglog(h, y_diff_central, 'r');
hold off;


