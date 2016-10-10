# -*- mode: python -*-

import os
import sys

block_cipher = None

regrid_dir = '../'
anaconda_dir = os.path.join(os.environ['HOME'], 'anaconda2')

# Edit this to specify the location to libllvmlite.so
libllvmlite = os.path.join(anaconda_dir, 'lib/python2.7/site-packages/llvmlite/binding/libllvmlite.so')
if not os.path.exists(libllvmlite):
    print("ERROR: can't find libllvmlite.so, please edit libllvm var in regrid.spec.")
    sys.exit(1)

a = Analysis([os.path.join(regrid_dir, 'regrid.py')],
             pathex=[regrid_dir],
             binaries=[(libllvmlite, 'llvmlite/binding')],
             datas=[],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='regrid',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='regrid')
