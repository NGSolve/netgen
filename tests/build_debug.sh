cd
mkdir -p build/netgen
cd build/netgen
cmake ../../src/netgen -DUSE_CCACHE=ON -DBUILD_TYPE=DEBUG
make -j12
make install
