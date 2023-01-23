function figure003_02_03(Y, p)

colors = {'r';'g';'b';'c';'m';'k'};

subplot(4,4,[1:4])
x = Y{1,2};
plot(x,Y{1,1},'k--')
hold on
for j = 2:size(Y,1)
    if j <= length(colors)+1
        c = colors{j-1};
    else
        c = colors{end};
    end
    for k = 1:size(Y{j,2},1)
        plot(Y{j,2}(k,:),Y{j,1}(k,:),'Color',c)
    end
end
hold off
axis tight
m_ylim = ylim;

subplot(4,4,[5:8])
[px,fx] = periodogram(Y{1,1} - mean(Y{1,1}),[],[],1);
px = filter(ones(1,7)/7, 1, px);
plot(fx,px);

ind = (9:16);
for j = 1:length(ind)
    if j >= size(Y,1)
        continue
    end
    if j < length(colors)
        c = colors{j};
    else
        c = colors{end};
    end
    subplot(4,4,ind(j))
    plot([])
    plot(Y{j+1,1}','Color',c)
    ylim(m_ylim);
    title([num2str(round(p(j)*100, 0)) '%'])
end

end