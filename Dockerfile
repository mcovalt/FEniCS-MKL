FROM ubuntu:16.04
SHELL ["/bin/bash", "-c"]
ARG intelkey
ARG uid
ARG gid
# Make sure args are set
RUN [ -n "$intelkey" ]
RUN [ -n "$uid" ]
RUN [ -n "$gid" ]
# Add user
RUN groupadd -r lol && useradd -r -g lol lol
# Generate locales
RUN apt-get -q update && apt-get -q install -y --no-install-recommends \
        locales \
 && rm -rf /var/lib/apt/lists/* \
 && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
 && locale-gen
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8
#=======================================================================================
# Versions
#=======================================================================================
##########
# FEniCS #
##########
ENV MSHR_VERSION="2017.2.0" \
    DOLFIN_VERSION="2017.2.0.post0" \
    DOLFINADJOINT_VERSION="2017.2.0" \
    FENICS_VERSION=">=2017.2.0,<2018.1.0"
#########
# Intel #
#########
ENV INTEL_VERSION_MAJOR="2018" \
    INTEL_VERSION_MINOR="3"
###########
# Compute #
###########
ENV PETSC_VERSION="3.8.3" \
    SLEPC_VERSION="3.8.2" \
    PETSC4PY_VERSION="3.8.1" \
    SLEPC4PY_VERSION="3.8.0" \
    IPOPT_VERSION="3.12.9" \
    NUMPY_VERSION="1.14.5" \
    SCIPY_VERSION="1.1.0"
################
# Supplemental #
################
ENV HDF5_VERSION="1.8.21" \
    SZIP_VERSION="2.1.1" \
    ZLIB_VERSION="1.2.11" \
    SWIG_VERSION="3.0.12" \
    MPI4PY_VERSION="3.0.0" \
    BOOST_VERSION="1.60.0" \
    PYBIND11_VERSION="2.2.1" \
    CYTHON_VERSION="0.28.2" \
    BOOST_VERSION="1.60.0" \
    PCRE_VERSION="8.42" \
    SYMPY_VERSION="1.1.1"
#=======================================================================================
# apt packages
#=======================================================================================
RUN apt-get -q update && apt-get -q install -y --no-install-recommends \
        bison \
        cmake \
        cpio \
        flex \
        g++ \
        gcc \
        git-core \
        libbz2-dev \
        libeigen3-dev \
        libfreetype6-dev \
        libgmp-dev \
        libmpfr-dev \
        libpng12-dev \
        libstdc++6 \
        nano \
        pkg-config \
        python \
        python3-dev \
        python3-pip \
        python3-setuptools \
        sudo \
        wget \
        zsh \
 && rm -rf /var/lib/apt/lists/*
#=======================================================================================
# Compile
#=======================================================================================
WORKDIR /tmp
############################
# Intel Parallel Studio XE #
############################
COPY intel.cfg /tmp/intel.cfg
RUN name="intel" \
 && url="http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/12998/parallel_studio_xe_${INTEL_VERSION_MAJOR}_update${INTEL_VERSION_MINOR}_cluster_edition_online.tgz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && echo "ACTIVATION_SERIAL_NUMBER=${intelkey}" >> /tmp/intel.cfg \
 && ./install.sh -s /tmp/intel.cfg \
    # Remove docs
 && rm -rf /opt/intel/documentation_${INTEL_MAJOR_VERSION} \
 && rm -rf /opt/intel/samples_${INTEL_MAJOR_VERSION} \
 && rm -rf /opt/intel/mkl/examples \
 && rm -rf /opt/intel/mkl/benchmarks \
 && rm -rf /opt/intel/compilers_and_libraries/linux/mpi/benchmarks \
    # Remove static MKL libs
 && rm -f /opt/intel/mkl/lib/intel64/libmkl_*.a \
 && rm -rf /tmp/* /var/tmp/*
###################
# Build variables #
###################
RUN echo "#!/bin/bash" > /usr/local/bin/vars.sh \
 && echo "export BUILD_THREADS=$( nproc )" >> /usr/local/bin/vars.sh \
 && echo "export INTEL_MAJOR_VERSION=2018" >> /usr/local/bin/vars.sh \
 && echo "source /opt/intel/parallel_studio_xe_\${INTEL_MAJOR_VERSION}/bin/psxevars.sh > /dev/null" >> /usr/local/bin/vars.sh \
 && echo "export CC=icc" >> /usr/local/bin/vars.sh \
 && echo "export CXX=icpc" >> /usr/local/bin/vars.sh \
 && echo "export FC=ifort" >> /usr/local/bin/vars.sh \
 && echo "export F77=ifort" >> /usr/local/bin/vars.sh \
 && echo "export F90=ifort" >> /usr/local/bin/vars.sh \
 && echo "export MPICC=mpiicc" >> /usr/local/bin/vars.sh \
 && echo "export MPICXX=mpiicpc" >> /usr/local/bin/vars.sh \
 && echo "export MPIF77=mpiifort" >> /usr/local/bin/vars.sh \
 && echo "export MPIF90=mpiifort" >> /usr/local/bin/vars.sh \
 && echo "export CFLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export CXXFLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export FCFLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export F77FLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export F90FLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export MPICCFLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export MPICXXFLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export MPIF77FLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export MPIF90FLAGS='-O3 -xHost -ip'" >> /usr/local/bin/vars.sh \
 && echo "export MKL_NUM_THREADS=1" >> /usr/local/bin/vars.sh \
 && echo "export OMP_NUM_THREADS=1" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
#########
# PETSc #
#########
RUN name="petsc" \
 && url="http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && ./configure \
        --with-cc=mpiicc \
        --with-cxx=mpiicpc \
        --with-fc=mpiifort \
        --with-shared-libraries \
        --with-debugging=0 \
        --with-c-support \
        --with-cxx-dialect=C++11 \
        --known-mpi-shared-libraries=1 \
        --COPTFLAGS="-O3 -xHost -ip" \
        --CXXOPTFLAGS="-O3 -xHost -ip" \
        --FOPTFLAGS="-O3 -xHost -ip" \
        --with-blas-include=/opt/intel/mkl/include \
        --with-blas-lib="-L/opt/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl" \
        --with-lapack-include=/opt/intel/mkl/include \
        --with-lapack-lib="-L/opt/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl" \
        --with-blacs-include=/opt/intel/mkl/include \
        --with-blacs-lib="-L/opt/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl" \
        --with-scalapack-include=/opt/intel/mkl/include \
        --with-scalapack-lib="-L/opt/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl" \
        --with-mkl_pardiso-dir=/opt/intel/mkl/ \
        --with-mkl_cpardiso-dir=/opt/intel/mkl/ \
        --download-mumps \
        --download-hypre \
        --download-metis \
        --download-parmetis \
        --download-suitesparse \
        --download-ptscotch \
        --download-superlu \
        --download-superlu_dist \
        --prefix=/opt/petsc \
 && make PETSC_DIR=/tmp/petsc PETSC_ARCH=arch-linux2-c-opt all \
 && make PETSC_DIR=/tmp/petsc PETSC_ARCH=arch-linux2-c-opt install \
 && echo "export PETSC_DIR=/opt/petsc" >> /usr/local/bin/vars.sh \
 && echo "export PETSC_ARCH=\"\"" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
#########
# SLEPc #
#########
RUN name="slepc" \
 && url="http://slepc.upv.es/download/distrib/slepc-${SLEPC_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && ./configure \
       --prefix=/opt/slepc \
 && make SLEPC_DIR=/tmp/slepc MAKE_NP=${BUILD_THREADS} \
 && make SLEPC_DIR=/tmp/slepc install \
 && echo "export SLEPC_DIR=/opt/slepc" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
########
# zlib #
########
RUN name="zlib" \
 && url="https://www.zlib.net/zlib-${ZLIB_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && ./configure \
        --prefix=/opt/zlib \
 && make -j ${BUILD_THREADS} && make install \
 && echo "export ZLIB_DIR=/opt/zlib" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
########
# szip #
########
RUN name="szip" \
 && url="https://support.hdfgroup.org/ftp/lib-external/szip/${SZIP_VERSION}/src/szip-${SZIP_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && ./configure \
        --prefix=/opt/szip \
 && make -j ${BUILD_THREADS} && make install \
 && echo "export SZIP_DIR=/opt/szip" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
########
# HDF5 #
########
RUN name="hdf5" \
 && url="https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION:0:3}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && CC=mpiicc ./configure \
        --enable-parallel \
        --with-zlib=${ZLIB_DIR} \
        --with-szlib=${SZIP_DIR} \
        --enable-build-mode=production \
        --prefix=/opt/hdf5 \
 && make -j ${BUILD_THREADS} && make install \
 && echo "export HDF5_DIR=/opt/hdf5" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
########
# PCRE #
########
RUN name="pcre" \
 && url="ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-${PCRE_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && ./configure \
        --disable-dependency-tracking \
        --enable-utf8 \
        --enable-pcre8 \
        --enable-pcre16 \
        --enable-pcre32 \
        --enable-unicode-properties \
        --enable-jit \
        --prefix=/opt/pcre \
 && make -j ${BUILD_THREADS} && make install \
 && echo "export PCRE_DIR=/opt/pcre" >> /usr/local/bin/vars.sh \
 && echo "export LD_LIBRARY_PATH=/opt/pcre/lib:\$LD_LIBRARY_PATH" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
########
# SWIG #
########
RUN name="swig" \
 && url="http://downloads.sourceforge.net/swig/swig-${SWIG_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && ./configure \
        --without-octave \
        --without-java \
        --with-pcre-prefix=${PCRE_DIR} \
        --prefix=/opt/swig \
 && make -j ${BUILD_THREADS} && make install \
 && echo "export SWIG_DIR=/opt/swig" >> /usr/local/bin/vars.sh \
 && echo "export PATH=/opt/swig/bin:\${PATH}" \
 && rm -rf /tmp/* /var/tmp/*
############
# pybind11 #
############
RUN name="pybind11" \
 && url="https://github.com/pybind/pybind11/archive/v${PYBIND11_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && mkdir build \
 && cd build \
 && cmake \
        -DPYBIND11_TEST=False \
        -DCMAKE_INSTALL_PREFIX=/opt/pybind11 \
        ../ \
 && make -j ${BUILD_THREADS} && make install \
 && echo "export PYBIND11_DIR=/opt/pybind11" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
##########
# Cython #
##########
RUN name="cython" \
 && url="https://pypi.python.org/packages/source/C/Cython/Cython-${CYTHON_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && LDSHARED='icc -shared' python3 setup.py install \
 && rm -rf /tmp/* /var/tmp/*
#########
# Numpy #
#########
RUN name="numpy" \
 && url="https://github.com/numpy/numpy/releases/download/v${NUMPY_VERSION}/numpy-${NUMPY_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && echo "[mkl]" >> site.cfg \
 && echo "library_dirs = /opt/intel/mkl/lib/intel64" >> site.cfg \
 && echo "include_dirs = /opt/intel/mkl/include" >> site.cfg \
 && echo "mkl_libs = mkl_rt" >> site.cfg \
 && echo "lapack_libs =" >> site.cfg \
 && sed -i "s/icc -m64 -fPIC -fp-model strict -O3/icc -m64 -fPIC -fp-model strict -O3 -xHost/g" numpy/distutils/intelccompiler.py \
 && sed -i "s/-fp-model strict -O1/-xhost -fp-model strict -fPIC/g" numpy/distutils/fcompiler/intel.py \
 && python3 setup.py \
        build \
            -j ${BUILD_THREADS} \
        config \
            --compiler=intelem build_clib \
            --compiler=intelem build_ext \
            --compiler=intelem \
        install \
 && rm -rf /tmp/* /var/tmp/*
#########
# Scipy #
#########
RUN name="scipy" \
 && url="https://github.com/scipy/scipy/releases/download/v${SCIPY_VERSION}/scipy-${SCIPY_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && python3 setup.py \
        config \
            --compiler=intelem \
            --fcompiler=intelem build_clib \
            --compiler=intelem \
            --fcompiler=intelem build_ext \
            --compiler=intelem \
            --fcompiler=intelem \
        install \
 && rm -rf /tmp/* /var/tmp/*
#########
# Boost #
#########
RUN name="../opt/boost" \
 && url="https://cytranet.dl.sourceforge.net/project/boost/boost/${BOOST_VERSION}/boost_$( echo ${BOOST_VERSION} | tr . _ ).tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && cd /opt/boost/tools/build \
 && ./bootstrap.sh \
        --with-toolset=intel-linux \
 && ./b2 \
        install \
        --prefix=/opt/boost \
 && export PATH=/opt/boost/bin:${PATH} \
 && cd /opt/boost \
 && b2 \
        -j ${BUILD_THREADS} \
        --with-filesystem \
        --with-iostreams \
        --with-math \
        --with-program_options \
        --with-system \
        --with-thread \
        --with-timer \
        --with-regex \
        --build-dir=/tmp/boost-build \
        threading=multi \
        variant=release \
        address-model=64 \
        toolset=intel \
        stage \
 && echo "export BOOST_DIR=/opt/boost" >> /usr/local/bin/vars.sh \
 && echo "export BOOST_ROOT=/opt/boost" >> /usr/local/bin/vars.sh \
 && echo "export PATH=/opt/boost/bin:\${PATH}" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
#########################
# other Python packages #
#########################
RUN source /usr/local/bin/vars.sh \
 && LDSHARED='icc -shared' pip3 --no-cache-dir install \
        matplotlib \
        --no-binary :all: \
 && pip3 --no-cache-dir install \
        fenics${FENICS_VERSION} \
        mpi4py==${MPI4PY_VERSION} \
        petsc4py==${PETSC4PY_VERSION} \
        ply \
        slepc4py==${SLEPC4PY_VERSION} \
        sympy==${SYMPY_VERSION} \
        --no-binary :all:
##########
# DOLFIN #
##########
RUN name="dolfin" \
 && url="https://bitbucket.org/fenics-project/dolfin/downloads/dolfin-${DOLFIN_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && mkdir build \
 && cd build \
 && CC=MPIICC cmake \
        -DCMAKE_BUILD_TYPE="Release" \
        -DCMAKE_CXX_COMPILER=mpiicpc \
        -DCMAKE_CXX_FLAGS="-O3 -xHost -ip -mkl -std=c++11" \
        -DCMAKE_C_COMPILER=mpiicc \
        -DCMAKE_C_FLAGS="-O3 -xHost -ip -mkl" \
        -DCMAKE_Fortran_COMPILER=mpiifort \
        -DCMAKE_Fortran_FLAGS="-O3 -xHost -ip -mkl" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DHDF5_C_COMPILER_EXECUTABLE=${HDF5_DIR}/bin/h5pcc \
        -DHDF5_hdf5_LIBRARY_RELEASE=${HDF5_DIR}/lib/libhdf5.so \
        -DHDF5_z_LIBRARY_RELEASE=${ZLIB_DIR}/lib/libz.so \
        -DPETSC_DIR=${PETSC_DIR} \
        -DSLEPC_DIR=${SLEPC_DIR} \
        -DDOLFIN_ENABLE_SPHINX=OFF \
        -DDOLFIN_ENABLE_UMFPACK=ON \
        -DDOLFIN_ENABLE_CHOLMOD=ON \
        -DZLIB_LIBRARY=${ZLIB_DIR}/lib/libz.a \
        -DZLIB_INCLUDE_DIR=${ZLIB_DIR}/include \
        -DDOLFIN_ENABLE_PASTIX=OFF \
        -DDOLFIN_ENABLE_TRILINOS=OFF \
        -DDOLFIN_ENABLE_BENCHMARKS=ON \
        -DDOLFIN_ENABLE_VTK=OFF \
        -DDOLFIN_ENABLE_QT=OFF \
        -DDOLFIN_AUTO_DETECT_MPI=ON \
        -DMPI_CXX_COMPILER=$( which mpiicpc ) \
        -DMPI_C_COMPILER=$( which mpiicc ) \
        -DMPI_Fortran_COMPILER=$( which mpiifort ) \
        -DHDF5_LIBRARIES=${HDF5_DIR}/lib/libhdf5.so \
        -DHDF5_C_COMPILER_EXECUTABLE=${HDF5_DIR}/bin/h5pcc \
        -DHDF5_hdf5_LIBRARY_RELEASE=${HDF5_DIR}/lib/libhdf5.so \
        -DPYTHON_EXECUTABLE:FILEPATH=$( which python3 ) \
        -DBOOST_ROOT=${BOOST_DIR} \
        -DSWIG_DIR=${SWIG_DIR} \
        -DSWIG_EXECUTABLE=${SWIG_DIR}/bin/swig \
        ../ \
 && make -j ${BUILD_THREADS} && make install \
 && echo "source /usr/local/share/dolfin/dolfin.conf" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
#########
# meshr #
#########
RUN name="meshr" \
 && url="https://bitbucket.org/fenics-project/mshr/downloads/mshr-${MSHR_VERSION}.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && mkdir build \
 && cd build \
 && cmake ../ \
 && make -j ${BUILD_THREADS} && make install \
 && rm -rf /tmp/* /var/tmp/*
#########
# Ipopt #
#########
RUN name="ipopt" \
 && url="https://www.coin-or.org/download/source/Ipopt/Ipopt-${IPOPT_VERSION}.tgz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && ./configure \
        --enable-shared \
        --with-blas-incdir=/opt/intel/mkl/include \
        --with-blas-lib="-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl" \
        --with-lapack-incdir=/opt/intel/mkl/include \
        --with-lapack-lib="-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl" \
        --prefix=/opt/ipopt \
 && make -j ${BUILD_THREADS} && make install \
 && echo "export IPOPT_DIR=/opt/ipopt" >> /usr/local/bin/vars.sh \
 && echo "export LD_LIBRARY_PATH=/opt/ipopt/lib:\$LD_LIBRARY_PATH" >> /usr/local/bin/vars.sh \
 && rm -rf /tmp/* /var/tmp/*
###########
# pyipopt #
###########
RUN name="pyipopt" \
 && url="https://github.com/xuy/pyipopt/archive/master.tar.gz" \
 && mkdir /tmp/$name && wget -nv -O - ${url} | tar xz -C /tmp/$name --strip-components=1 && cd /tmp/$name \
 && source /usr/local/bin/vars.sh \
    # Build
 && sed -i "s#/usr/local#/opt/ipopt#g" setup.py \
 && sed -i "s/coinmumps/dmumps/g" setup.py \
 && sed -i "s/coinblas/mkl_rt/g" setup.py \
 && sed -i "s/coinlapack/mkl_rt/g" setup.py \
 && sed -i "s/coinmetis/metis/g" setup.py \
 && sed -i "s#library_dirs=\[IPOPT_LIB\]#library_dirs=[IPOPT_LIB,'/opt/petsc/lib']#g" setup.py \
 && LDSHARED='icc -shared' python3 setup.py build \
 && python3 setup.py install \
 && rm -rf /tmp/* /var/tmp/*
#############
# pyadjoint #
#############
RUN source /usr/local/bin/vars.sh \
 && pip3 --no-cache-dir install \
        git+https://bitbucket.org/dolfin-adjoint/pyadjoint.git@${DOLFINADJOINT_VERSION} \
        --no-binary :all:
##=======================================================================================
## User environment
##=======================================================================================
RUN echo "lol:lol" | chpasswd \
 && echo "lol ALL=(ALL) NOPASSWD: ALL" | sudo EDITOR='tee -a' visudo \
 && mkdir /home/lol \
 && mkdir /home/lol/shared \
 && usermod -u $uid lol \
 && groupmod -g $gid lol \
 && chown -R lol:lol /home/lol
VOLUME /home/lol/shared
WORKDIR /home/lol/shared
USER lol
RUN wget https://github.com/robbyrussell/oh-my-zsh/raw/master/tools/install.sh -O - | /usr/bin/zsh || true \
 && sed -i "s/robbyrussell/avit/g" /home/lol/.zshrc \
 && git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting \
 && sed -i "s/  git/  git zsh-syntax-highlighting/g" /home/lol/.zshrc \
 && echo "source /usr/local/bin/vars.sh" >> /home/lol/.zshrc
CMD /usr/bin/zsh