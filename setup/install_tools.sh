#!/usr/bin/env bash
# Install Java 17 + Nextflow 26.10.4 into ~/tools exactly as on the original machine (conf/env.sh expects this).
# Prerequisites (apt): docker.io (or Docker CE), curl, wget, python3, poppler-utils. Add yourself to the docker group.
set -euo pipefail
mkdir -p ~/tools && cd ~/tools
curl -fsSL -o jdk17.tar.gz "https://api.adoptium.net/v3/binary/latest/17/ga/linux/x64/jdk/hotspot/normal/eclipse"
mkdir -p jdk17 && tar -xzf jdk17.tar.gz -C jdk17 --strip-components=1 && rm jdk17.tar.gz
mkdir -p nextflow_bin && cd nextflow_bin
curl -fsSL https://get.nextflow.io -o nextflow_installer
JAVA_HOME=~/tools/jdk17 PATH=~/tools/jdk17/bin:$PATH NXF_HOME=~/tools/.nextflow NXF_VER=26.10.4 bash nextflow_installer
chmod +x nextflow && rm -f nextflow_installer
JAVA_HOME=~/tools/jdk17 NXF_HOME=~/tools/.nextflow NXF_VER=26.10.4 ./nextflow -version | grep -i version
