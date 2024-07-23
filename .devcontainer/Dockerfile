FROM mcr.microsoft.com/devcontainers/universal:2

RUN apt update && apt upgrade
RUN apt install fish && chsh fish


ENV PYODIDE_VERSION=0.26.1
ENV EMSDK_VERSION=3.1.58
ENV PYTHON_VERSION=3.12.1

# get emscripten
WORKDIR /root
RUN git clone https://github.com/emscripten-core/emsdk.git
RUN cd emsdk && ./emsdk install ${EMSDK_VERSION} && ./emsdk activate ${EMSDK_VERSION}

RUN pip install pyodide-build==${PYODIDE_VERSION} --break-system-packages

RUN pyodide venv env26