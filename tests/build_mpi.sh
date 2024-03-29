cd
mkdir -p build/netgen
cd build/netgen
cmake ../../src/netgen \
  -DCHECK_RANGE=ON \
  -DUSE_CCACHE=ON \
  -DUSE_MPI=ON \
  -DUSE_OCC=OFF \
  -DENABLE_UNIT_TESTS=ON
make -j12
make install
