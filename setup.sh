DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH=$DIR/python:$PYTHONPATH
export LD_LIBRARY_PATH=$DIR/analyzer/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$DIR/analyzer/lib:$DYLD_LIBRARY_PATH
export PATH=$DIR/bin:$PATH
export HEPMC2ROOT_PATH=$DIR
