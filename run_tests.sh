#!/bin/bash

PYTHONPATH=../asp/:$PYTHONPATH

echo PYTHONPATH
echo ${PYTHONPATH}

if [ -z "${PYTHON}" ]
then
    PYTHON=python
fi
if [ -z "${PYTHONARGS}" ]
then
    PYTHONARGS=
fi

PYTHONPATH=`pwd`:${PYTHONPATH} ${PYTHON} ${PYTHONARGS} tests/1d_fft_test.py

PYTHONPATH=`pwd`:${PYTHONPATH} ${PYTHON} ${PYTHONARGS} tests/fft_test.py
