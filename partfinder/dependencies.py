# Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# PartitionFinder also includes the PhyML program, the RAxML program, and the
# PyParsing library, all of which are protected by their own licenses and
# conditions, using PartitionFinder implies that you agree with those licences
# and conditions as well.


import imp

DEPENDENCIES = ['numpy', 'pandas', 'tables', 'pyparsing', 'scipy', 'sklearn']

for dep in DEPENDENCIES:
    try:
        imp.find_module(dep)
    except ImportError:
        print('\n\n\n **** ERROR **** \n')
        print('Could not find the dependency %s, please check that you have '
	    	  'followed the installation instructions in the manual '
	    	  ' and try again.' %(dep))
        print('\nPartitionFinder 2 (unlike PartitionFinder 1) requires a few'
        	  ' other Python packages to work. All of these '
			  'can be installed very easily, in a single click, by installing '
			  'the Python 2.7.x '
			  'version of the Anaconda Python distrubition, which you can '
			  'find here (remember to get the 2.7 version!): '
			  '\n\nhttps://www.continuum.io/downloads ')
        print('\nFurther instructions are in the installation section of the '
			  'manual, which points out a number of other differences between '
			  'PartitionFinder 1 and PartitionFinder 2 too.')
        print('\n *************** \n\n\n')

        raise ImportError