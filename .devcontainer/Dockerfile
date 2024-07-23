FROM archlinux

RUN pacman -Syu --noconfirm
RUN pacman -S git gcc cmake make openblas fish python nodejs python-pip --noconfirm


ENV PYODIDE_VERSION=0.26.1
ENV EMSDK_VERSION=3.1.58
ENV PYTHON_VERSION=3.12.1

# get emscripten
WORKDIR /root
RUN git clone https://github.com/emscripten-core/emsdk.git
RUN cd emsdk && ./emsdk install ${EMSDK_VERSION} && ./emsdk activate ${EMSDK_VERSION}

RUN pip install pyodide_build==0.26.0 --break-system-packages

# RUN pyodide venv env26
