# FEniCS MKL
A Dockerfile for building an optimized python3 FEniCS stack using the Intel compilers and the Intel MKL. This includes FEniCS, meshr, and dolfin-adjoint with support for minimization via Ipopt and SciPy.

## Why use this Dockerfile?
1. It is a low effort way to speed up your existing FEniCS projects.
2. It is an insight into what is required when compiling your own FEniCS stack with the Intel compilers.

## Build Instructions
Three arguments must be defined at build time:
* `intelkey`: a serial number used to register the copy of Intel Parallel Studio XE within the Docker image. Intel provides this for free to students, classroom educators, and open-source contributors. See [here](https://software.intel.com/en-us/parallel-studio-xe/choose-download) for more details.
* `uid`: the user ID of the user who wishes to use this image. Most often this is stored in the environment variable `UID`.
* `gid`: the group ID of the user who wishes to use this image. Most often this is stored in the environment variable `GID`.

Assuming the UID and GID environment variables are set, enter the directory where the Dockerfile and `intel.cfg` are located and execute the following.

```bash
$ docker build . -t fenics-mkl --build-arg intelkey="XXXX-XXXXXXXX" --build-arg uid="$UID" --build-arg gid="$GID"
```

Expect it to take a while. There's a lot of code to download and compile.

## Usage Instructions

Initialize the Docker container by running the image and mapping a chosen directory for the image to access.

```bash
$ docker run -v '/a/local/directory':'/home/lol/shared' -it fenics-intel:latest
```

Exit the image then find the container ID.

```bash
$ docker ps -a
```

Then whenever you'd like to enter the container, start and attach to it. For example, if your container ID is `df6652b339cf`, execute the following.

```bash
$ docker start df6652b339cf && docker attach df6652b339cf
```

It is not required to use the container, but since the FFC cache for your FEniCS project is stored in the container, reusing the container will allow the JIT compiler to use previously compiled code.