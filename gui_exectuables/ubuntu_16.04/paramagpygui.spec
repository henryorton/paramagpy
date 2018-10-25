# -*- mode: python -*-

block_cipher = None


a = Analysis(['../../paramagpy/gui.py'],
             pathex=['/home/u5376227/Dropbox/PhD/git/paramagpy/gui_exectuables/ubuntu_16.04'],
             binaries=[],
             datas=[],
             hiddenimports=["numpy.lib.recfunctions"],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='paramagpygui',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=True )
