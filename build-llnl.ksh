#!/bin/ksh -f

case "$SYS_TYPE" in
    toss_3_x86_64_ib)  for=toss3 ;;
    toss_3_aarch64_ib) for=aarch64 ;;
    blueos*)           for=pwr9 ;;
    *)                 for=default ;;
esac

while [ "$#" != 0 ]; do
    case "$1" in
        -clean)
            make clean
            exit 0 ;;
        -for)
            if [ -r "defs/config.for.$2" ]; then
                for="$2"
            else
                echo build: no such file as defs/config.for.$2
                exit 1
            fi
            shift ;;
        -help)
            echo "Usage:  build [options] [package1 [package2 ...]]"
            echo "options:"
            echo "          -for what    use file defs/config.for.what to configure" ;;
        *)
            echo "unknown argument" ;;
    esac
    shift
done

. "./defs/config.for.$for"

export CC CFLAGS LDFLAGS LIB

# Always make clean because the code builds so quickly. 
# Note that a "make clean" is required whenever switching architectures
# e.g. x86_64 to Blue Gene/Q
make clean

make fftpf3d

