FROM debian:bookworm-slim AS shimmer

RUN apt update

RUN apt install -y apt-utils build-essential gdb cmake cmake-curses-gui \
                   git vim libeigen3-dev liblua5.4-0 liblua5.4-dev \
                   libsqlite3-dev sqlite3 sqlite3-tools libboost-all-dev

RUN useradd --create-home shimmer
USER shimmer

WORKDIR /home/shimmer/
RUN git clone --recursive https://github.com/shimmerhydrogen/shimmer.git
RUN mkdir -p /home/shimmer/shimmer/shimmer++/build/
WORKDIR /home/shimmer/shimmer/shimmer++/build/
RUN ls -la ..
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc)

