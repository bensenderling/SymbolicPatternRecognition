function y = create003_01_01
% This code generates theoretical data to test the sumbolic pattern
% recognition algorithm.

% Number of different segments.
N = 4;
% Number of each segment.
M = 3;
% Data length of each segment.
n = 2000;
% Time step.
f = 100;
% Time variable for a single period.
t = (0:5*n-1)'/f;

s = [];
for i = 1:M
    s = [s;ceil(N*rand(5,1))];
end
    
y = [];
for i = 1:length(s)
    switch s(i)
        case 1
            y = [y;sin(2*pi*t(1:n))];
        case 2
            y = [y;sin(2*3*pi*t(1:n))];
        case 3
            y = [y;square(2*pi*t(1:n))];
        case 4
            y = [y;sawtooth(2*pi*t(1:n))];
    end
end

y = y + randn(length(y),1)*0;
