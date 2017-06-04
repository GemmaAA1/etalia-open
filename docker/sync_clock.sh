#!/usr/bin/env bash

# Run this to fix the "box out of sync" bug that causes communication with S3 to fail
# (fix found here: https://forums.docker.com/t/time-in-container-is-out-of-sync/16566/8)
docker run -it --rm --privileged --pid=host etalia/dev-ssh nsenter -t 1 -m -u -n -i date -u $(date -u +%m%d%H%M%Y)
docker run -it --rm --privileged --pid=host etalia/dev nsenter -t 1 -m -u -n -i date -u $(date -u +%m%d%H%M%Y)
