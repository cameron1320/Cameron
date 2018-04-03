function orbitAnimation(obj)

distance = .5*obj.Fs/100;
threshold = .5*(max(obj.refchop(:,ceil(end/2))) - mean(obj.refchop(:,ceil(end/2))));
[~, loc] = findpeaks(obj.ref,'MinPeakHeight',threshold,'MinPeakDistance',distance);
pksx=obj.X(loc);
pksy=obj.Y(loc);
Fsh = obj.Fs;
maxlm = max(obj.X);
speedarray = (obj.Speed(end)-obj.Speed(1))/length(pksy)*(0:length(pksy)-1)+obj.Speed(1);
for i = 1:1:length(pksx)-1;
    plot(obj.X(loc(i):loc(i+1)-fix((1/10)*(loc(i+1)-loc(i)))),obj.Y(loc(i):loc(i+1)-fix((1/10)*(loc(i+1)-loc(i)))),'-k',pksx(i),pksy(i),'.','markersize',15)
    text = [' Running Speed: ',num2str(fix(speedarray(i)),4),' (RPM)'];
    uicontrol('Style', 'text',...
        'String',text,...
        'Units','normalized',...
        'Position', [.25 .9 0.5 .05]);
    axis([-maxlm maxlm -maxlm maxlm]);
    axis square
    xlabel('Horizontal Amplitude (mils)');
    ylabel('Vertical Amplitude (mils)');
%     pause(1/Fsh)
    pause

end
end