function ax = fig_raster(data,ops,t,cell_oi,ind_front,ind_rear)


	% count spikes time course
	ops.posterior_t_edges = (-2:0.2:6)+t;
	[count,ops] = classifier.count_spk.time_course(data,ops);
	Posterior = ops.classifier.predict(ops.Mdl,count',ops);

	% plot
	ax = np(3,1);

	% raster
	toi = [-1 5] + t;
	tmp = cellfun(@(x) x(find(x>toi(1) & x<toi(2))) , data.spikes(cell_oi),'UniformOutput',false);
	plt.raster_smooth(tmp(ind_front),toi,ax(1),'kernel_width',200); 
	plt.raster_smooth(tmp(ind_rear), toi,ax(2),'kernel_width',200); 

	% posterior
	Colors = [1 0 0;0 0 0];
	for jj = 1:2
		h = plot(ax(3),ops.posterior_t,Posterior(:,jj),'Color',Colors(jj,:),'LineWidth',0.7);
		area(ax(3),ops.posterior_t,Posterior(:,jj),'FaceColor',Colors(jj,:),'EdgeColor','none','FaceAlpha',0.1);
	end

	% figure setting
	% arrayfun(@(ii) set(ax(ii),), 1:2);
	arrayfun(@(ii) set(ax(ii),'XLim',[-1 5]+t,'Visible','off'), 1:3);
	colorbar off;
	set(gcf,'Position',[0 0 100 300]);
	ef;
	
end