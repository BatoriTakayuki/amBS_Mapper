FROM python:3.12-slim-bookworm
ENV DEBIAN_FRONTEND=noninteractive PYTHONUNBUFFERED=1 PYTHONPATH=/app AMBS_DATA=/data MPLCONFIGDIR=/tmp/matplotlib
RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools samtools tabix ca-certificates curl fonts-ipafont-gothic build-essential zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt \
    && apt-get purge -y build-essential && apt-get autoremove -y
COPY ambs_mapper /app/ambs_mapper
COPY tests /app/tests
WORKDIR /data
EXPOSE 8501
ENTRYPOINT ["python", "-m", "ambs_mapper"]
CMD ["web", "--host", "0.0.0.0"]
