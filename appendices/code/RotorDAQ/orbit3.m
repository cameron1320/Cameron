
function orbit3(obj)
N = length(obj.X(:));
F = (obj.Speed(end) - obj.Speed(1))/N*(0:N-1) + obj.Speed(1);
s = fix(sqrt(N));
Xlin = obj.X(1:s^2);
Ylin = obj.Y(1:s^2);
w = F(1:s^2);
x = reshape(Xlin,s,s);
y = reshape(Ylin,s,s);
o = reshape(w,s,s);
C = sqrt(x.^2+y.^2);
h = mesh(o,x,y,C);
axis([1000 2000 -inf inf -inf inf]);
xlabel('Running Speed (RPM)');
ylabel('Horizontal Amplitude (mils)');
zlabel('Vertical Amplitude (mils)');
h.FaceColor = 'none';
h.MeshStyle = 'column';
end