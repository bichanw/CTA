function fig_line_shade(ax,x,y,color)
plot(ax,x,y,'Color',color,'LineWidth',0.7);
area(ax,x,y,'FaceColor',color,'EdgeColor','none','FaceAlpha',0.1);
			