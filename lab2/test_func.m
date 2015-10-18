clear all,close all
keys=[]
for i = 1:8
path=[ '/Users/JackChen/Desktop/ASAS/lab2/test00' num2str(i) '.wav' ];

key=getPhoneNum(path);
keys=[keys '\n' key];

end


fprintf([keys '\n'] )
%%
load TA_answer_partial.mat
keys=answer_map_partial.keys;
% fileID = fopen('exp.txt','w');
% keys(ii)
% answer_map_partial(keys{ii})

for ii = 1:length(keys)
    fprintf([keys{ii} ': ' answer_map_partial(keys{ii}) '\n']);
end
% fclose(fileID);