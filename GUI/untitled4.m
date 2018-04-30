for i =0:1:11
    ReadData('250418-2-Matth-2/dataEEG2504Rest' num2str(i+1) '.txt',3,125,5,i*10);
end