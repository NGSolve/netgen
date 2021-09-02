cd
mkdir -p build/netgen
cd build/netgen
cmake \
    -DUSE_CCACHE=ON \
    -DBUILD_TYPE=DEBUG \
    -DENABLE_UNIT_TESTS=ON \
    -DUSE_OCC=ON \
    -DCHECK_RANGE=ON \
    -DUSE_CGNS=ON \
    -DCMAKE_INSTALL_PREFIX=/usr \
    ../../src/netgen
make -j12
make install
