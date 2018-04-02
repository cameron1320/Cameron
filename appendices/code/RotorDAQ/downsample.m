function downsample(obj, r)
obj.X = downsample(obj.X, r);
obj.Y = downsample(obj.Y, r);
obj.ref = downsample(obj.ref, r);
obj.Fs = obj.Fs/r;
end