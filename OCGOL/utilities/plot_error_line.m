function [h] = plot_error_line(range_mean_sem,style_l,width_l,color_l)
%input range, mean, and sem at 3xtime matrix
%input style, width, and color
%return handle to the plot

h = errorbar(range_mean_sem(1,:),range_mean_sem(2,:),range_mean_sem(3,:),'LineStyle',style_l,'Linewidth',width_l,'Color',color_l);


end

