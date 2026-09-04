FROM debian:bookworm-slim AS builder

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        libhts-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build/isoumi
COPY Makefile VERSION ./
COPY src ./src

RUN make -C src \
    && strip src/isoumi

FROM python:3.12-slim-bookworm

ARG ISOUMI_VERSION=0.1.1
LABEL org.opencontainers.image.title="IsoUMI" \
      org.opencontainers.image.description="Isoform-aware long-read UMI correction with JupyterLab" \
      org.opencontainers.image.version="${ISOUMI_VERSION}" \
      org.opencontainers.image.source="https://github.com/dawangran/IsoUMI" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        ca-certificates \
        libgomp1 \
        samtools \
        tini \
    && rm -rf /var/lib/apt/lists/* \
    && python -m pip install --no-cache-dir "jupyterlab>=4.2,<5"

COPY --from=builder /build/isoumi/src/isoumi /usr/local/bin/isoumi
COPY docker/start-jupyter.sh /usr/local/bin/start-isoumi-jupyter
COPY examples/minimal.sam /opt/isoumi/examples/minimal.sam
COPY notebooks /workspace/notebooks

RUN chmod 0755 /usr/local/bin/isoumi /usr/local/bin/start-isoumi-jupyter \
    && mkdir -p /workspace

ENV JUPYTER_ROOT_DIR=/workspace \
    JUPYTER_PORT=8888 \
    PYTHONUNBUFFERED=1

WORKDIR /workspace
EXPOSE 8888

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["/usr/local/bin/start-isoumi-jupyter"]
