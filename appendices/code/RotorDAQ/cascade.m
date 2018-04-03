function cascade(obj)
figure('name','Cascade Plot')
waterfall(obj.freq, obj.Speed, obj.Zf)
xlabel('Frequency of Vibration (Hz)');
ylabel('Running Speed (RPM)');
zlabel('Amplitude of Vibration (mils)');
axis([-100 100 -inf inf 0 inf])
end