# femtoscopy

```bash
git clone https://github.com/lbavinh/femtoscopy.git
cd femtoscopy
mkdir build
cd build
cmake ..
make -j4
./PicoDstFemtoscopy -i  <input ROOT file or list> -o Pi.root -cme 3
root -l -b -q ../macros/PlotCF.cpp
```