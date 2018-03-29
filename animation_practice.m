x=0:0.05:4*pi;
y=sin(x);
curve=animatedline;
set(gca,'XLim',[0 4*pi],'YLim',[-1 1]);
for i=1:length(x)
    addpoints(curve,x(i),y(i));
    drawnow
end