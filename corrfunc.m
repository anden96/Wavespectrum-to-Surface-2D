function corrFunc = corrfunc(surface)

z = surface
n = length(z);
delta = zeros(length(n),1);

for i = 2:1:length(delta)
    delta(i) = rms(n(i-1) - n(i));
    
    end
    deltalog = 10*log10(delta);
    deltalogsort = sort(deltalog);
    
    plot(linspace(0,1,length(deltalogsort)),deltalogsort)
    
    corrFunc = deltalogsort;
    
end