function [cmap_return] = return_paper_colormap()
%returns the normalized RGB values of the colormap used for the episodic
%memory/splitter paper
%A - royal blue
%B - crimson red
%A&B - dark magenta
%A&B - A - light blue (light sky blue)
%A&B - B - light red (light coral red)

cmap_return = [65,105,225;
                220,20,60;
                139, 0, 139;
                135,206,250;
                240,128,128]./255;


end

