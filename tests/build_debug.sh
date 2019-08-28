cd
mkdir -p build/netgen
cd build/netgen
cmake ../../src/netgen -DUSE_CCACHE=ON -DBUILD_TYPE=DEBUG -DENABLE_UNIT_TESTS=ON -DUSE_OCC=ON
make -j12
make install
