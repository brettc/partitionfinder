import sys
print sys.path

from distutils.core import setup
import py2app

setup(
    app=['pfgui.py'],
    data_files=[
    	'Resources/blob.tiff',
    	'Resources/blob_file.icns',
    	'Resources/blob_app.icns'
    ],
    options = dict(
    	py2app = dict(
    		#includes = ['objc', 'Foundation', 'AppKit'],
    		#packages = ['GUI'],
    		plist = dict(
    			#CFBundleSignature = "BLBE", # If you are using mac_creator
    			CFBundleIconFile = "blob_app",
    			CFBundleDocumentTypes = [
    				dict(CFBundleTypeName = "PartitionFinder Document",
    					CFBundleTypeRole = "Editor",
    					#CFBundleTypeOSTypes = ["BLOB"], # If you are using mac_type
    					CFBundleTypeExtensions = ["cfg"],
    					CFBundleTypeIconFile = "blob_file",
    					)]
    			),
    	),
    )
)
