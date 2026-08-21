#!/bin/bash
# R08 transient VM. Same pattern as the 1.25 TiB pull: scratch disk, ship logs to GCS every
# 60s, run the job, write a sentinel. Teardown is guaranteed externally by maxRunDuration +
# instanceTerminationAction=DELETE, because a minimal self-delete role was refused earlier.
set -uo pipefail
BUCKET=gs://heydonto-quantara-lungcdx
RUN=$BUCKET/data-request-2026-08/_run/r08
LOG=/var/log/r08.log
exec > >(tee -a "$LOG") 2>&1
echo "=== R08 VM STARTUP $(date -u) ==="

mkdir -p /mnt/scratch
if [ -b /dev/nvme0n1 ]; then
  mkfs.ext4 -F /dev/nvme0n1 >/dev/null 2>&1 && mount /dev/nvme0n1 /mnt/scratch
fi
chmod 1777 /mnt/scratch
df -h /mnt/scratch

export DEBIAN_FRONTEND=noninteractive
command -v python3 >/dev/null || { apt-get update -y >/dev/null; apt-get install -y python3 >/dev/null; }
python3 -c "import numpy" 2>/dev/null || {
  apt-get update -y >/dev/null
  apt-get install -y python3-numpy >/dev/null || pip3 install --break-system-packages numpy
}
python3 --version; python3 -c "import numpy;print('numpy',numpy.__version__)"

( while true; do gcloud storage cp "$LOG" "$RUN/log.txt" >/dev/null 2>&1; sleep 60; done ) &
SHIP=$!

gcloud storage cp "$RUN/vm_job_r08.py" /opt/vm_job_r08.py
WORKERS=12 python3 /opt/vm_job_r08.py
RC=$?
echo "=== job rc=$RC $(date -u) ==="

kill $SHIP 2>/dev/null
gcloud storage cp "$LOG" "$RUN/log.txt" >/dev/null 2>&1
echo "$RC" | gcloud storage cp - "$RUN/DONE_rc_${RC}" >/dev/null 2>&1
# best-effort self-delete; the Google-enforced cap is the real guarantee
NAME=$(curl -s -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/name)
ZONE=$(curl -s -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/zone | awk -F/ '{print $NF}')
gcloud compute instances delete "$NAME" --zone "$ZONE" --quiet || echo "self-delete refused; relying on maxRunDuration"
