FROM cyversevice/jupyterlab-scipy:latest

USER root

# Install R and related packages.
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list \
        && apt-key adv --keyserver keys.gnupg.net --recv-key '51716619E084DAB9' \
        && apt-get update \
        && apt-get install -y r-base-dev libcairo2-dev libxt-dev xvfb \
