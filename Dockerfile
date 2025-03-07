FROM mambaorg/micromamba:1.4.9

# metadata labels
LABEL base.image="mambaorg/micromamba:1.4.9"
LABEL dockerfile.version="1"
LABEL software="vibecheck"
LABEL software.version="2025.01.07"
LABEL description="Conda environment for Vibecheck. Vibecheck: Software package for assigning O1 Vibrio cholerae genome sequences to canonical lineages."
LABEL website="https://github.com/watronfire/Vibecheck"
LABEL license="GNU General Public License v3.0"
LABEL license.url="https://github.com/watronfire/Vibecheck/blob/main/LICENSE"
LABEL maintainer="Nathaniel L. Matteson"
LABEL maintainer.email="nmatteson@bwh.harvard.edu"

WORKDIR /home/mambauser/
COPY --chown=$MAMBA_USER:$MAMBA_USER . ./vibecheck/
WORKDIR /home/mambauser/vibecheck/

USER root
RUN apt-get update && apt-get install -y --no-install-recommends git
USER mambauser
RUN micromamba install -y -n base -f environment.yaml \
    && micromamba clean --all --yes
# Specify path to gain access to conda environment without terminal model
ENV PATH="/opt/conda/envs/base/bin:/opt/conda/bin:/opt/conda/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

RUN pip install .
RUN vibecheck -v

WORKDIR /data/
