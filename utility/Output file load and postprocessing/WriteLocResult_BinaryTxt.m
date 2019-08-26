function WriteLocResult_BinaryTxt(WriteDat, FileName)

% WriteDat = LocArry;

% TotalFrame = LocArry(end,end);

% FileName='LocArry_whole.txt';

% FileName = sprintf('LocArry_whole_%df.txt', TotalFrame);

WriteDat = WriteDat';
WriteDat = single(WriteDat);

fid = fopen(FileName,'w');
fwrite(fid,WriteDat(:),'single');

fclose(fid);

