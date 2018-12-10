cd
mkdir -p build/netgen
cd build/netgen
cmake ../../src/netgen -DUSE_CCACHE=ON -DENABLE_CPP_CORE_GUIDELINES_CHECK=ON
make -j12
make install
