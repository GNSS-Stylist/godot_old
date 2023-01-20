scons --clean -j10 platform=windows target=release_debug
scons --clean -j10 platform=windows target=release_debug float=64
scons --clean -j10 platform=windows tools=no target=release_debug float=64
scons --clean -j10 platform=windows tools=no target=release float=64

scons -j10 platform=windows target=release_debug
scons -j10 platform=windows target=release_debug float=64
scons -j10 platform=windows tools=no target=release_debug float=64
scons -j10 platform=windows tools=no target=release float=64

pause
