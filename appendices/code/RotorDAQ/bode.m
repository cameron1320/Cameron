function bode(obj)
figure('Name','Bode Plot')
subplot(211)
plot(obj.Speed,obj.Amp.XAmp,obj.Speed,obj.Amp.YAmp)
xlabel('Running Speed (RPM)');
ylabel('Amplitude (mils)');
legend('Horizontal Plane','Vertical Plane');
subplot(212)
plot(obj.Speed,obj.Phase.PhaseX,obj.Speed,obj.Phase.PhaseY)
xlabel('Running Speed (RPM)');
ylabel('Phase lag (deg.)');
set(gca,'ytick',-720:45:720,'YTickLabel',{'0','45','90','135','180','225','270','315'})
end