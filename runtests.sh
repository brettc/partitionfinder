# Not much here.
# But it is easier to type "sh runtests.sh" than this...
# --buffer captures output unless a test fails
# tests is the folder where our tests are
# python -m unittest discover --buffer --verbose tests
nosetests -v tests
