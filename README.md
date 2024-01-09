# ASC-ODE
A package for solving ordinary differential equations


TODO:
* add your basic-linalg as git submodule
* adjust CMakeLists.txt to include your marix.h / vector.h files
* adjust matrix/vector interface

## Running mass_spring

```sh
python3.12 -m venv venv
source venv/bin/activate

pip install jupyterlab numpy numpy pybind11
pip install pythreejs
pip install BLA/. -v

# Now build the ASC-ODE cmake project with VS Code or with something like this:
# cmake -S. -B build
# cd build && make build && cd ..

mv build/mass_spring/*.so build/mass_spring/mass_spring.so
jupyter-lab
```