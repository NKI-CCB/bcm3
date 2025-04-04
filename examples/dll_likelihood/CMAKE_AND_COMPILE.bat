mkdir build
cd build
cmake ..
msbuild likelihood.sln /p:Configuration=Release
cd ..
