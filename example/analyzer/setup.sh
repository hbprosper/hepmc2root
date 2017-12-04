
if [ "$PYTHONPATH" == "" ]; then
    export PYTHONPATH=$PWD/python
else
    export PYTHONPATH=$PWD/python:$PYTHONPATH
fi

if [ "$LD_LIBRARY_PATH" == "" ]; then
    export LD_LIBRARY_PATH=$PWD/lib
else
    export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
fi
