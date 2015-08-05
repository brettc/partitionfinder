# Get the submodules
git submodule init
git submodule update

# Building Phyml
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

mv -f ./submodules/phyml/src/phyml ./programs

