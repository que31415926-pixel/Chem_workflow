# Cloud Build receives a user-authorized Linux x86-64 ORCA 6 OpenMPI 4.1.6 installer.
FROM debian:bookworm-slim AS licensed
RUN apt-get update && apt-get install -y --no-install-recommends xz-utils bzip2 && rm -rf /var/lib/apt/lists/*
COPY orca.run /tmp/installer.run
RUN chmod 0700 /tmp/installer.run && /tmp/installer.run -- -p /opt/orca && test -x /opt/orca/orca

FROM mambaorg/micromamba:2.0.5
USER root
RUN micromamba install -y -n base -c conda-forge python=3.12 openmpi=4.1.6 crest=3.0.2 "xtb>=6.7,<7" rdkit "numpy>=1.26,<3" "scipy>=1.12,<2" python-libsbml=5.20.5 google-cloud-storage && micromamba clean --all --yes
COPY --from=licensed /opt/orca /opt/orca
WORKDIR /opt/chemflow
COPY chemflow2 ./chemflow2
COPY vendor ./vendor
ENV PATH=/opt/orca:/opt/conda/bin:${PATH} \
    LD_LIBRARY_PATH=/opt/orca:/opt/conda/lib \
    PYTHONPATH=/opt/chemflow \
    CHEMFLOW_ORCA=/opt/orca/orca \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
RUN test "$(uname -m)" = x86_64 && test -x /opt/orca/orca_scf_mpi && \
    ldd /opt/orca/orca_scf_mpi > /tmp/orca-libraries && ! grep -q "not found" /tmp/orca-libraries && \
    mpirun --version | grep -F 4.1.6 && mpirun --oversubscribe -np 2 hostname && \
    crest --version && xtb --version && python3 -c "import scipy, rdkit, libsbml, google.cloud.storage, chemflow2"
