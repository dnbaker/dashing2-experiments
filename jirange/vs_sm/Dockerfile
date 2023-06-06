FROM fedora:32

RUN sudo dnf upgrade --refresh -y \
 && sudo dnf install gcc-c++ cmake procps-ng which python3-pip wget xz git zlib-devel unzip findutils -y

RUN pip3 install --upgrade pip

# Get Dashing 2
RUN wget -O /tmp/dashing2_s128.xz -q https://github.com/dnbaker/dashing2-binaries/raw/main/linux/v2.1.14/dashing2_s128.xz
RUN xz -d /tmp/dashing2_s128.xz
RUN sudo mv /tmp/dashing2_s128 /bin/dashing2 && sudo chmod a+x /bin/dashing2
RUN dashing2 --version || true

# Get Dashing 1
RUN wget -O /tmp/dashing_s128.xz -q https://github.com/dnbaker/dashing-binaries/raw/main/linux/dashing_s128.xz
RUN xz -d /tmp/dashing_s128.xz
RUN sudo mv /tmp/dashing_s128 /bin/dashing && sudo chmod a+x /bin/dashing
RUN dashing --version || true

# Get Sourmash
RUN pip3 install sourmash
RUN sourmash --version

# Get Mash
RUN wget -O /tmp/mash-Linux64-v2.3.tar -q https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar
RUN cd /tmp && tar xvf mash-Linux64-v2.3.tar && sudo mv mash-Linux64-v2.3/mash /bin/mash && sudo chmod a+x /bin/mash
RUN mash --version

# Get BinDash
RUN cd /tmp && git clone https://github.com/zhaoxiaofei/bindash.git \
 && mkdir -p bindash/release && cd bindash/release \
 && cmake -DCMAKE_BUILD_TYPE=Release .. && make \
 && sudo mv bindash /bin/bindash && sudo chmod a+x /bin/bindash

# Get FastANI
RUN wget https://github.com/ParBLiSS/FastANI/releases/download/v1.33/fastANI-Linux64-v1.33.zip \
 && unzip fastANI-Linux64-v1.33.zip \
 && sudo mv fastANI /bin/fastANI && chmod a+x /bin/fastANI \
 && rm -f fastANI-Linux64-v1.33.zip

# Clean up
RUN rm -rf /tmp/*
