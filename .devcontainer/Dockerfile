FROM ubuntu:22.04

# -------------- Base C++ development -------------

RUN apt update
RUN apt install -y git curl zip openssh-client
RUN	apt install -y gcc g++ cmake

# packages required by c++ packages installed through vcpkg
RUN apt install -y pkg-config bison python3 python3-distutils libxmu-dev libxi-dev libgl-dev libxt-dev libglvnd-dev liburing-dev
RUN apt-get clean && apt autoremove -y

# -------------- VCPKG -----------------------

RUN git clone https://github.com/microsoft/vcpkg.git /opt/vcpkg
RUN /opt/vcpkg/bootstrap-vcpkg.sh
ENV VCPKG_ROOT=/opt/vcpkg
ENV VCPKG_TOOLCHAIN=/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
ENV VCPKG_FEATURE_FLAGS=versions
ENV PATH=$PATH:/opt/vcpkg