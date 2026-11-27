# Cloud Build receives the user-authorized ORCA 6.1.0 Linux x86-64 installer linked against Open MPI 4.1.8
# (exact file and SHA-256 in toolchain.lock.json; cloud_setup refuses any other installer). conda-forge's 4.x series ends at 4.1.6, so
# OpenMPI 4.1.8 is built from the official source release and placed first on PATH/LD_LIBRARY_PATH.
FROM debian:bookworm-slim AS licensed
RUN apt-get update && apt-get install -y --no-install-recommends xz-utils bzip2 && rm -rf /var/lib/apt/lists/*
COPY orca.run /tmp/installer.run
RUN chmod 0700 /tmp/installer.run && /tmp/installer.run -- -p /opt/orca && test -x /opt/orca/orca

FROM debian:bookworm-slim AS mpi
ARG OMPI_VERSION=4.1.8
RUN apt-get update && apt-get install -y --no-install-recommends build-essential gfortran curl ca-certificates bzip2 \
    && rm -rf /var/lib/apt/lists/*
RUN curl -fsSL "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-${OMPI_VERSION}.tar.bz2" -o /tmp/ompi.tar.bz2 \
    && tar -xjf /tmp/ompi.tar.bz2 -C /tmp && cd /tmp/openmpi-${OMPI_VERSION} \
    && ./configure --prefix=/opt/openmpi --disable-mpi-fortran --without-verbs > /dev/null \
    && make -j"$(nproc)" > /dev/null && make install > /dev/null && rm -rf /tmp/openmpi-* /tmp/ompi.tar.bz2

FROM mambaorg/micromamba:2.0.5
USER root
# Open MPI (built from source above) needs an rsh/ssh agent binary to exist even for single-node launches;
# without it mpirun aborts with "plm_rsh_agent: ssh : rsh ... could not be found" (build b1230b1b, 2026-11-26).
RUN apt-get update && apt-get install -y --no-install-recommends openssh-client && rm -rf /var/lib/apt/lists/* && command -v ssh
# Exact versions = toolchain.lock.json (the versions resolved by build b1230b1b on 2026-11-26).
# Change them only together with toolchain.lock.json; v2/tests/test_toolchain_lock.py checks they agree.
RUN micromamba install -y -n base -c conda-forge python=3.12.16 crest=3.0.2 xtb=6.8.0 tblite=0.8.0 rdkit=2026.09.2 numpy=2.5.5 scipy=1.18.1 python-libsbml=5.20.5 google-cloud-storage=3.17.2 && micromamba clean --all --yes
COPY --from=mpi /opt/openmpi /opt/openmpi
COPY --from=licensed /opt/orca /opt/orca
WORKDIR /opt/chemflow
COPY chemflow2 ./chemflow2
COPY vendor ./vendor
ENV PATH=/opt/orca:/opt/openmpi/bin:/opt/conda/bin:${PATH} \
    LD_LIBRARY_PATH=/opt/orca:/opt/openmpi/lib:/opt/conda/lib \
    PYTHONPATH=/opt/chemflow \
    CHEMFLOW_ORCA=/opt/orca/orca \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_plm=isolated
# OMPI_MCA_plm=isolated: every ORCA job runs on one VM, so mpirun forks ranks locally and never needs ssh/rsh
# (verified on the Mac with Open MPI 4.1.6 by the DeepSeek diagnosis of build b1230b1b). openssh-client above is a fallback.
# Self-check (fix 2026-11-26). ORCA 6.1 no longer ships a separate orca_scf_mpi module, so the check no
# longer names one module. Every check runs even if an earlier one fails, and all results are printed,
# so a single build reports every problem at once. The build fails if any check fails.
#  - every ORCA *_mpi module resolves all its libraries, and the MPI modules link /opt/openmpi
#  - mpirun is 4.1.8 and can start 2 ranks
#  - CREST, xTB and the Python stack import
#  - a real 2-process ORCA HF/def2-SVP water job terminates normally
RUN fail=0; \
    check() { if sh -c "$2" > /tmp/selfcheck.log 2>&1; then echo "SELFCHECK PASS $1"; else echo "SELFCHECK FAIL $1"; tail -n 25 /tmp/selfcheck.log | sed 's/^/    | /'; fail=1; fi; }; \
    check arch 'test "$(uname -m)" = x86_64'; \
    check orca-exe 'test -x /opt/orca/orca && ls -l /opt/orca/orca'; \
    check orca-mpi-modules 'ls /opt/orca/ | grep -E "_mpi$" | tee /dev/stderr | grep -q .'; \
    check orca-mpi-libraries 'rc=0; for b in /opt/orca/*_mpi; do if ldd "$b" | grep -q "not found"; then echo "$b:"; ldd "$b" | grep "not found"; rc=1; fi; done; exit $rc'; \
    check orca-links-openmpi418 'n=0; for b in /opt/orca/*_mpi; do ldd "$b" | grep -q "/opt/openmpi/lib/libmpi" && n=$((n+1)); done; echo "$n MPI modules link /opt/openmpi"; test "$n" -gt 0'; \
    check mpirun-version 'mpirun --version | tee /dev/stderr | grep -F 4.1.8'; \
    check mpirun-2-ranks 'mpirun --oversubscribe -np 2 hostname'; \
    check crest 'crest --version'; \
    check xtb 'xtb --version'; \
    check python-stack 'python3 -c "import scipy, rdkit, libsbml, google.cloud.storage, chemflow2; print(\"python ok\")"'; \
    check orca-parallel-water 'mkdir -p /tmp/orcatest && cd /tmp/orcatest && printf "! HF def2-SVP\n%%pal nprocs 2 end\n* xyz 0 1\nO 0 0 0.117\nH 0 0.757 -0.467\nH 0 -0.757 -0.467\n*\n" > w.inp && /opt/orca/orca w.inp > w.out 2>&1; tail -n 15 w.out; grep -q "ORCA TERMINATED NORMALLY" w.out'; \
    grep -i -E "parallel|MPI-process" /tmp/orcatest/w.out 2>/dev/null | head -n 3 | sed 's/^/SELFCHECK INFO /'; \
    rm -rf /tmp/orcatest /tmp/selfcheck.log; \
    if [ "$fail" -ne 0 ]; then echo "SELFCHECK RESULT: FAILED (see FAIL lines above)"; exit 1; fi; \
    echo "SELFCHECK RESULT: ALL PASSED"
