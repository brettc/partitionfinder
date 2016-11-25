
clean-test:
	find ./tests -name analysis -and -d | xargs rm -rf 
	find ./tests -name log.txt | xargs rm

test:
	py.test --verbose --color=yes tests

quick-test:
	py.test -v -k "not full" -x

init-repo:
	# Get the submodules (raxml and phyml)
	git submodule init
	git submodule update

# This needs work
# Building Phyml ----------------------------------
# This is not functioning
# phyml: init-repo ./programs/phyml
# 	cd submodules/phyml
# 	aclocal
# 	autoheader
# 	autoconf -f
# 	automake -f --add-missing
# 	./configure
# 	make clean
# 	make
# 	cd ../..
# 	mv -f ./submodules/phyml/src/phyml ./programs

# Building Raxml ----------------------------------
# We use the SSE as it is the most reliable
# raxml: init-repo ./programs/raxml
# 	cd submodules/raxml
# 	make -f Makefile.SSE3.gcc
# 	cd ../..
# 	mv -f ./submodules/raxml/raxmlHPC-SSE3 ./programs/raxml

.PHONY: init-repo clean-test test quick-test
