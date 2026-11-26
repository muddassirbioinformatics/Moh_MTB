# Our own image for summary parsing, step PDFs and the clinical report. Rebuild with:
#   docker build -t mtb-report:1.1.0 containers/mtb-report
# 1.1.0: + pypdf (merged per-sample step PDF)
FROM python:3.12.7-slim-bookworm
RUN pip install --no-cache-dir reportlab==4.2.5 pyyaml==6.0.2 matplotlib==3.9.2 pypdf==5.1.0 \
 && apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
ENV MPLCONFIGDIR=/tmp/mpl
