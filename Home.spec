# -*- mode: python ; coding: utf-8 -*-

a = Analysis(
    ['run.py'],
    pathex=[],
    binaries=[],
    datas=[
        ('assets/', 'assets/'),
        ('example-data/', 'example-data/'),
        ('.streamlit/', '.streamlit/'),
    ],
    hiddenimports=[
        'streamlit',
        'streamlit.runtime',
        'streamlit.web.cli',
    ],
    hookspath=['./hooks'],
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='Corromics',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
)

# Bundle for macOS
app = BUNDLE(
    exe,
    name='Corromics.app',
    icon='assets/corromics_icon.ico',
    bundle_identifier='com.functionalmetabolomics.corromics',
)
