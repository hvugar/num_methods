setenv('MW_MINGW64_LOC','E:\MinGW\7.3.0\x86_64');

cd E:\project\hvugar\num_methods\trunk\optimal\problem2H\matlab\

if libisloaded('problem2H')
    unloadlibrary problem2H
end

copyfile ..\..\bin\problem2H.dll .
copyfile ..\..\bin\minimum.dll .
copyfile ..\..\problem2H\problem2h_exporter.h .
copyfile ..\..\problem2H\problem2h_global.h .

loadlibrary('problem2H', 'problem2h_exporter.h');

calllib('problem2H', 'init_pr');
calllib('problem2H', 'setPenaltyR', 1.0);

disp('Initialized.');