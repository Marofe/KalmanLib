function draw_gaussian(mu,P,color,width)
t=-4*P:0.01:4*P;
f=1/(sqrt(2*pi)*P)*exp(-(mu-t).^2/(2*P^2));
plot(t,f,color,'linewidth',width)
end