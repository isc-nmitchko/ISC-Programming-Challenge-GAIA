ARG IMAGE=intersystems/iris-community:latest-em
FROM $IMAGE

WORKDIR /home/irisowner/dev
COPY . .

## Embedded Python environment
ENV IRISUSERNAME="_SYSTEM"
ENV IRISPASSWORD="SYS"
ENV IRISNAMESPACE="USER"
ENV PYTHON_PATH=/usr/irissys/bin/
ENV PATH="/usr/irissys/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/irisowner/bin"

## Install processing libraries into IRIS embedded Python and expose the src modules
ENV PYTHONPATH=/home/irisowner/dev/src
RUN pip3 install --no-cache-dir --target /usr/irissys/mgr/python polars

## Pre-build the OPT-IN libdeflate C scan kernel used by gaia_fast.py (portable -O3,
## no -march=native so it stays safe across CPUs). Done as root because COPY'd files
## are root-owned. The default gaia.py is pure polars and does NOT need this; it only
## benefits the gaia_fast.py benchmarking path, which also lazily recompiles the .so
## at runtime (volume-mount case) and falls back to polars if gcc/libdeflate are gone.
USER root
RUN gcc -O3 -fopenmp -shared -fPIC src/gaiascan.c -o src/gaiascan.so \
        /lib/x86_64-linux-gnu/libdeflate.so.0 -lm \
    && echo "built gaiascan.so" || echo "gaiascan.so build skipped; gaia_fast.py uses polars fallback"
USER irisowner

RUN --mount=type=bind,src=.,dst=. \
    iris start IRIS && \
	iris merge IRIS merge.cpf && \
	iris session IRIS < iris.script && \
    iris stop IRIS quietly safely
