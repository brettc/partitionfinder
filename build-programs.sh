# Get the submodules (raxml and phyml)
git submodule init
git submodule update

# Building Phyml ----------------------------------
# This is modified from the confphy...
# It assumes that you've done:
# brew install automake autoconf libtool 
cd submodules/phyml
aclocal
autoheader
autoconf -f
# Note: brew renames this...
glibtoolize
automake -f --add-missing
./configure
make clean
make

cd ../..
mv -f ./submodules/phyml/src/phyml ./programs


# Building Raxml ----------------------------------
cd submodules/raxml
make -f submodules/raxml/Makefile.SSE3.mac

cd ../..
mv -f ./submodules/raxml/raxmlHPC-SSE3 ./programs/raxml
