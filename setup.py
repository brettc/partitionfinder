import sys
print sys.path

from distutils.core import setup
import py2app

setup(
    app=['wxgui.py'],
    # app=['PartitionFinder.py'],
    # data_files=[
        # 'Resources/blob.tiff',
        # 'Resources/blob_file.icns',
        # 'Resources/blob_app.icns'
    # ],
    options = dict(
    	py2app = dict(
    		#includes = ['objc', 'Foundation', 'AppKit'],
    		#packages = ['GUI'],
            # resources=['resources/License.txt'],
            iconfile='resources/pf.icns',
            packages='wx',
            site_packages=True,
    		plist = dict(
                CFBundleName               = "PartitionFinder",
                CFBundleShortVersionString = "1.0",     # must be in X.X.X format
                CFBundleGetInfoString      = "PartitionFinder 1.0",
                CFBundleExecutable         = "PartitionFinder",
                CFBundleIdentifier         = "com.robertlanfear.partitionfinder",
                # CFBundleIconFile = "blob_app",
                # CFBundleDocumentTypes = [
                    # dict(CFBundleTypeName = "PartitionFinder Document",
                        # CFBundleTypeRole = "Editor",
                        # #CFBundleTypeOSTypes = ["BLOB"], # If you are using mac_type
                        # CFBundleTypeExtensions = ["cfg"],
                        # CFBundleTypeIconFile = "blob_file",
                        # )]
            ),
        ),
    )
)

# import py2app

    # Build the .app file
    # setup(
        # options=dict(
            # py2app=dict(
                # iconfile='resources/myapp-icon.icns',
                # packages='wx',
                # site_packages=True,
                # resources=['resources/License.txt'],
                # plist=dict(
                    # CFBundleName               = "MyAppName",
                    # CFBundleShortVersionString = "0.2.5",     # must be in X.X.X format
                    # CFBundleGetInfoString      = "MyAppName 0.2.5",
                    # CFBundleExecutable         = "MyAppName",
                    # CFBundleIdentifier         = "com.example.myappname",
                # ),
            # ),
        # ),
        # app=[ 'myapp.py' ]
