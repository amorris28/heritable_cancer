#!/bin/bash

REMOTE="nird"
REMOTE_DIR="~/heritable_cancer"
LOCAL_RESULTS_DIR="./"

echo "==> Syncing files to server..."
rsync -av --checksum --times \
    tsd_run.sh \
    run_pipeline.py \
    main.nf \
    nextflow.config \
    $REMOTE:$REMOTE_DIR/

rsync -av --checksum --times src/ $REMOTE:$REMOTE_DIR/src/

if [[ -f hazrd_latest.tar.gz ]]; then
    echo "==> Syncing hazrd_latest.tar.gz to server..."
    rsync -av --checksum hazrd_latest.tar.gz $REMOTE:$REMOTE_DIR/hazrd_latest.tar.gz
else
    echo "==> No hazrd_latest.tar.gz found locally — skipping (copy it here manually to update hazrd on TSD)"
fi

echo "==> Submitting job p2521 (blocking until complete)..."
ssh -t $REMOTE "cd $REMOTE_DIR && python run_pipeline.py --project p2521"

echo "==> Submitting job p212 (blocking until complete)..."
ssh -t $REMOTE "cd $REMOTE_DIR && python run_pipeline.py --project p212"

echo "==> Copying results to local..."
scp $REMOTE:$REMOTE_DIR/public_p2521.zip $LOCAL_RESULTS_DIR
scp $REMOTE:$REMOTE_DIR/public_p212.zip $LOCAL_RESULTS_DIR
unzip -o public_p2521.zip
unzip -o public_p212.zip

echo "==> Done! Results are in $LOCAL_RESULTS_DIR"
