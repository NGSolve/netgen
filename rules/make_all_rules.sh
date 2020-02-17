g++ makerlsfile.cpp -o makerls
./makerls hexa.rls        ../libsrc/meshing/hexarls.cpp       hexrules
./makerls prisms2.rls     ../libsrc/meshing/prism2rls.cpp     prismrules2
./makerls pyramids.rls    ../libsrc/meshing/pyramidrls.cpp    pyramidrules
./makerls pyramids2.rls   ../libsrc/meshing/pyramid2rls.cpp   pyramidrules2
./makerls quad.rls        ../libsrc/meshing/quadrls.cpp       quadrules
./makerls tetra.rls       ../libsrc/meshing/tetrarls.cpp      tetrules
./makerls triangle.rls    ../libsrc/meshing/triarls.cpp       triarules
rm makerls

# node: prisms2rls is currently not used (prism2rls_2.cpp is not compiled)
# ./makerls prisms2.rls     ../libsrc/meshing/prism2rls_2.cpp   prismrules2
