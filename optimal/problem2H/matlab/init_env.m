setenv('MW_MINGW64_LOC','E:\MinGW\5.3.0_x86_64');

cd E:\project\hvugar\num_methods\trunk\optimal\problem2H\matlab\

if libisloaded('problem2H')
    unloadlibrary problem2H
end

copyfile ..\..\bin\problem2H.dll .
copyfile ..\..\bin\minimum.dll .
copyfile ..\..\problem2H\exporter.h .

loadlibrary('problem2H', 'exporter.h');

calllib('problem2H', 'init_pr');
calllib('problem2H', 'setPenaltyR', 1.0);

disp('Initialized.');