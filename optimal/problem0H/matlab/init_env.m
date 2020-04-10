disp('Starting init_env...');
setenv('MW_MINGW64_LOC','E:\MinGW\7.3.0\x86_64');

cd E:\project\hvugar\num_methods\trunk\optimal\problem0H\matlab\

if libisloaded('problem0H')
    unloadlibrary problem0H
end

copyfile ..\..\bin\minimum.dll .
copyfile ..\..\bin\imaging.dll .
copyfile ..\..\bin\problem0H.dll .
copyfile ..\..\problem0H\problem0h_exporter.h .
copyfile ..\..\problem0H\problem0h_global.h .

loadlibrary('problem0H', 'problem0h_exporter.h');

calllib('problem0H', 'init_problem');
calllib('problem0H', 'setPenaltyR', 1.0);

disp('Initialized.');