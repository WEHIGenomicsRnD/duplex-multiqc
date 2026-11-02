FROM community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc

# Release tag of the plugin to install (e.g. v0.1.0)
ARG VERSION

LABEL org.opencontainers.image.source="https://github.com/WEHIGenomicsRnD/duplex-multiqc"
LABEL org.opencontainers.image.version="${VERSION}"

# Install from the GitHub release tarball (no git needed in the base image)
RUN test -n "${VERSION}" \
    && pip install --no-cache-dir \
        "https://github.com/WEHIGenomicsRnD/duplex-multiqc/archive/refs/tags/${VERSION}.tar.gz"
