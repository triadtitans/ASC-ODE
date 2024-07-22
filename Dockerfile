FROM ubuntu:24.04 AS pyodide_wheel

RUN apt-get update && apt-get -y upgrade
RUN apt-get update && apt-get -y install \
        nodejs \
        ccache \
        python3-pip \
        python3 \
        git \
        make \
        pkg-config \
        g++ \
        lbzip2 \
        xz-utils \
        autoconf \
        libtool \
        unzip \
        xxd \
        wget \
        python3.12-venv

ENV PYODIDE_VERSION=0.26.1
ENV EMSDK_VERSION=3.1.58
ENV PYTHON_VERSION=3.12.1

# get emscripten
WORKDIR /root
RUN git clone https://github.com/emscripten-core/emsdk.git
RUN cd emsdk && ./emsdk install ${EMSDK_VERSION} && ./emsdk activate ${EMSDK_VERSION}

RUN pip install pyodide-build==${PYODIDE_VERSION} --break-system-packages

RUN pyodide venv env26

WORKDIR /root/ASC-ODE
COPY . .

# time to build

SHELL ["/bin/bash", "-c"]
RUN source /root/emsdk/emsdk_env.sh && source /root/env26/bin/activate && pyodide build

WORKDIR /root/ASC-ODE/dist
RUN source /root/env26/bin/activate && \
        pip install lib_rigid_body-0.0.1-cp312-cp312-pyodide_2024_0_wasm32.whl --force-reinstall && \
        python3 -c "import lib_rigid_body as rb; print(dir(rb))"




FROM pyodide_wheel AS jupyterlite_build

WORKDIR /root/ASC-ODE
RUN pip install -r requirements.txt --break-system-packages

WORKDIR /root
RUN git clone https://github.com/triadtitans/rigid_body_interactive.git
WORKDIR /root/rigid_body_interactive

RUN TZ=UTC jupyter lite build --pyodide https://ngsolve.org/files/pyodide-0.26.0/ngsolve_pyodide.tar.bz2
RUN cp /root/ASC-ODE/dist/* /root/rigid_body_interactive/dist/static/pyodide


FROM jupyterlite_build AS server

EXPOSE 8000
WORKDIR /root/rigid_body_interactive/dist

