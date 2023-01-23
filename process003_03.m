function process003_03
% process003_02_03
% Inputs  - data loaded from a structure.
% Outputs - figures and a spreadsheet of the summary measures.
% Remarks
% - The code below will run the symbolic pattern recognition on joint
%   angles or synthetic data. The synthetic data is constructed from
%   randomly ordering sequences of sine, square and sawtooth waves.
% - The method uses a fourier spectrum to identify the dominant frequency
%   in the time series. It then uses this calculated a continuous wavelet
%   spectrum around that frequency. The location of the dominant frequency
%   is then used to find the most prominant sequence. This is done because
%   biological time series do not always have a constant frequency
%   spectrum.
% - Once a sequence is found all of its iterations are found. This is used
%   to find a mean sequence which is then used to again find all the
%   sequences. This continues untill the mean sequences from one iteration
%   to the other are minimally different.
% - Once a sequence is found the remaining data is sent into the pattern
%   recognition again. This process is done recursively untill too few data
%   remains to do further calculations.
% Future Work
% - The algorithm needs a stopping point. Right now this is 1000 data
%   points. At some point the data becomes too short to do a frequence or
%   wavelet spectrum on. This threshold likely varies for each time series.
%   Potential this limit could be based on the ability of the frequency and
%   wavelet spectrums to find a prominant sequence instead of a set data
%   length.
% Nov 2021 - Created by Ben Senderling, bensenderling@gmail.com
%%

% Use the synthetic data.
y = create003_01_01;

% Does the main work creating the histogram.
[Y, p] = sequenceHistogram(y);

% Create the figure for the histogram.
H = figure('visible','on');
figure003_02_03(Y, p)
drawnow

for i = 1:length(p)
    fprintf('%.0f of the data is accounted for in sequence one\n', p(i)*100)
end


end

function [Y, p] = sequenceHistogram(y)

y = (y-min(y))/range(y);

% Total number of symbols used.
N = 26;
% Starting character. 65 is A.
N2 = 65;
        % The time series is normalized from N2 to N2+N and rounded to whole
        % numbers so they can be converted to symbols.
        y2 = round(y*N+N2);
        % Convert the numbers to symbols.
        Y2 = char(y2);
        % This counts how many repetitions of a symbol are in a row. It cooresponds
        % to the compressed symbol sequence.
        Y2n = ones(length(Y2),1);
        % Creates a compressed symbol sequence.
        i = 2;
        while i < length(Y2)
            if Y2(i) == Y2(i-1)
                % Removed repeted characters.
                Y2(i) = [];
                % Increments value if a symbol is repeated.
                Y2n(i-1) = Y2n(i-1)+1;
                % Removed a spot for repeated characters.
                Y2n(end) = [];
            else
                i = i + 1;
            end
        end
        % Converts the symbolic sequence back into a numerical array. The speed was
        % increased by 1 s by removing this from the for loop.
        Y = {double(Y2)};
        Y{1,2} = (1:length(Y2))';
        
        error = 2;
        Y = sequenceFinder(Y, error);
        
        % Reorder the histogram from highest to lowest.
        [n,~] = cellfun(@size,Y(:,1));
        [~,I] = sort(n,'descend');
        Y = Y(I,:);
        
        % calculate the percentage of data accounted for.
        p = zeros(size(Y,1)-1,1);
        for i = 2:size(Y,1)
            p(i-1) = numel(Y{i,1})/numel(Y{1,1});
        end
        
end

function Y = sequenceFinder(Y, error)

% Select last sequence.
y = Y{end,1};
indexes = Y{end,2};
L = size(Y,1);

% Find prominant frequency to identify sequence length.
try
    [px,fx] = periodogram(y-mean(y),[],[],1);
    px = filter(ones(1,7)/7, 1, px);
    [~,I1] = max(px);
    sequenceLength = round(1/fx(I1));
catch
    Y = {};
    return
end

% Find a random start.
% start = ceil(rand(1)*(length(y) - sequenceLength));
% if start < 1
%     Y = {};
%     return
% end
% Find a start based on wavelet spectrum.
% The cwt() function has limits on it's frequency bounds. We want to use
% tight bounds to speed it up but not tighter than the limits.
f_a = 5;
f_diff = fx(I1)*(2^(0.1*f_a) - 1)/(1 + 2^(0.1*f_a));
f_lower = fx(I1) - f_diff;
f_upper = fx(I1) + f_diff;
% Calculate theoretical bounds and compare to calculated.
[minfreq,maxfreq] = cwtfreqbounds(length(y));
if f_lower < minfreq
    f_upper = f_upper + (minfreq - f_lower);
    f_lower = minfreq;
end
% Calculate the wavelet spectrum and find the index where the most
% prominant frequency has the most power.
[wt,f] = cwt(y-mean(y),'FrequencyLimits',[f_lower,f_upper]);
[~,I2] = min(abs(f - fx(I1)));
[~,start] = max(abs(wt(I2,ceil(sequenceLength/2) + 1:end - floor(sequenceLength/2) - 1)));
start = start + ceil(sequenceLength/2);

% Initial motif.
if rem(sequenceLength,2) == 0
    motif = y(start - sequenceLength/2:start + sequenceLength/2 - 1);
else
    motif = y(start - floor(sequenceLength/2):start + floor(sequenceLength/2));
end

m_diff = inf;
ind2 = 2;
while length(m_diff) == 1 || abs(m_diff(end) - m_diff(end-1)) > 1
    
    Y{L+1,1} = [];
    Y{L+1,2} = [];
    
    matches = zeros(length(y),1);
    
    % Find all sequences of the motif.
    k = 1;
    ind1 = 1;
    while k < length(y)-(sequenceLength)
        % Similar sequences are concidered within an error range.
        if sum(abs(y(k:k + sequenceLength - 1) - motif(:,ind2-1))) <= ceil(error*sequenceLength)
            Y{L+1,1}(ind1,:) = y(k:k + sequenceLength - 1);
            Y{L+1,2}(ind1,:) = indexes(k:k + sequenceLength - 1);
            ind1 = ind1 + 1;
            matches(k:k + sequenceLength - 1) = 1;
            k = k + sequenceLength;
        else
            k = k + 1;
        end
    end
    
    % If an matches were found find the difference with the previous
    % iteration and save the average sequence.
    if sum(matches) > 0
        m_diff(ind2,1) = sum(abs(mean(Y{2,1},1) - motif(:,ind2 - 1)'));
        motif(:,ind2) = mean(Y{2,1},1);
    else
        motif(:,ind2) = motif(:,ind2 - 1);
        m_diff(ind2,1) = m_diff(ind2-1,1);
    end
    ind2 = ind2 + 1;
    
end

% Find remaining data.
ind2 = find(matches == 0, 1);
% leftOver = {};
leftOver = {[],[]};
ind4 = 1;
while ind2 < length(matches) - 1
    ind3 = find(matches(ind2:end-1) == 0 & diff(matches(ind2:end)) == 1, 1) + ind2 - 1;
%     leftOver{ind4,1}(1,:) = y(ind2:ind3);
    leftOver{1} = [leftOver{1};y(ind2:ind3)];
%     leftOver{ind4,2}(1,:) = indexes(ind2:ind3);
    leftOver{2} = [leftOver{2};indexes(ind2:ind3)];
    ind4 = ind4 + 1;
    ind2 = find(matches(ind3 + 1:end) == 0, 1) + ind3;
end

for k = 1:size(leftOver,1)
    % Only rerun sequenceFinder if the segment is over a set length.
    if length(leftOver{k}) > 1000
        Y2 = sequenceFinder(leftOver(k,1:2), error);
        if ~isempty(Y2)
            Y(size(Y,1)+1:size(Y,1)+size(Y2,1)-1,1:2) = Y2(2:end,:);
        end
    end

end

end

