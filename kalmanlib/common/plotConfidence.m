function plot=plotConfidence(x,ciY,color)
plot=fill([x flip(x)],[ciY(1,:) flip(ciY(2,:))],color,'FaceAlpha', 0.2,'linestyle','none');
end