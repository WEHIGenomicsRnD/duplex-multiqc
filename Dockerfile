FROM community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc

LABEL org.opencontainers.image.source="https://github.com/WEHIGenomicsRnD/duplex-multiqc"

# The Wave image keeps its env under pixi and has no pip; put it on PATH for
# build steps and for Nextflow, which bypasses the ENTRYPOINT
ENV PATH="/opt/wave/.pixi/envs/default/bin:${PATH}"

# Install the plugin from the checkout (the release tag in CI). There is no git
# in the image, so the source is copied in. Dependencies (multiqc, matplotlib)
# already come with the base image.
COPY . /tmp/plugin
RUN python -m ensurepip \
    && python -m pip install --no-cache-dir --no-deps /tmp/plugin \
    && python -m pip uninstall -y pip \
    && rm -rf /tmp/plugin

# Fail the build if MultiQC or the plugin can't be loaded
RUN multiqc --version && python -c "import duplex_multiqc"
