# build and run

```
# clone this repository and its submodules
$ git clone --recurse-submodules https://github.com/w3ntao/smallpt-cpu.git

# build
$ cd smallpt-cpu
$ mkdir build; cd build
$ cmake ..; make -j

# render
$ ./smallpt-cpu
rendering (64 spp) took 10.600 seconds.
image saved to `smallpt_cpu_64.png`
```
