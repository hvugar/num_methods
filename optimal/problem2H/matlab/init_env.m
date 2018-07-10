%setenv('MW_MINGW64_LOC','E:\MinGW\5.3.0_x86_64');

cd E:\project\hvugar\num_methods\trunk\optimal\problem2H\matlab\

if libisloaded('problem2H')
    unloadlibrary problem2H
end

copyfile ..\..\bin\problem2H.dll .
copyfile ..\..\bin\minimum.dll .
copyfile ..\..\problem2H\exporter.h .

x0 = [+2.3400, -2.7400, +1.5800, +1.9500, +0.5000, -0.4000, -0.3000, +0.6000, +0.5500, +0.1400, +0.7400, +0.3700, +0.2800, +0.7500, +0.8500, +0.8900];

loadlibrary('problem2H', 'exporter.h');

calllib('problem2H', 'init_pr');
calllib('problem2H', 'setPenaltyR', 1);

disp('Initialized.');