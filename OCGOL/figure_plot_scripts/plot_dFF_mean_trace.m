function [s1,s2] = plot_dFF_mean_trace(traceA,traceB)
%plot A trials
sem_input =  @(x) nanstd(x,1)./sqrt(size(x,1));

s1 = shadedErrorBar(1:size(traceA,2),traceA  ,{@mean,sem_input},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s1.edge,'LineWidth',0.1,'LineStyle','-','Color',[[65,105,225]/255, 0.2]) %last value add transparency value
s1.mainLine.LineWidth = 2;
s1.mainLine.Color = [65,105,225]/255;
s1.patch.FaceColor = [65,105,225]/255;


%plot B trials
s2 = shadedErrorBar(1:size(traceB,2),traceB  ,{@mean,sem_input},'lineprops','-','transparent',true,'patchSaturation',0.20);
set(s2.edge,'LineWidth',0.1,'LineStyle','-','Color',[[ 220,20,60]/255, 0.2]) %last value add transparency value
s2.mainLine.LineWidth = 2;
s2.mainLine.Color = [ 220,20,60]/255;
s2.patch.FaceColor = [ 220,20,60]/255;
end

