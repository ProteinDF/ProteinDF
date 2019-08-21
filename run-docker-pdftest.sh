#!/bin/bash -eux

PDF_BUILDER_VER="latest"
#BRANCH=develop
PDF_BRIDGE_BRANCH="master"
PDF_PYTOOLS_BRANCH="master"
#DOCKER_TERM="-e COLUMNS=$COLUMNS -e LINES=$LINES -e TERM=$TERM"
DOCKER_CONTAINER_NAME="pdf-builder"

# -----------------------------------------------------------------------------
# test
# -----------------------------------------------------------------------------
run_test()
{
    docker exec -it \
       --env OMP_NUM_THREADS=4 \
       --env OMP_SCHEDULE=dynamic \
       ${DOCKER_CONTAINER_NAME} \
       pdf-check.sh --branch develop --workdir /tmp/pdf-check \
       serial 2>&1 | tee out.test_serial

    docker exec -it \
       --env OMP_NUM_THREADS=4 \
       --env OMP_SCHEDULE=dynamic \
       ${DOCKER_CONTAINER_NAME} \
       pdf-check.sh --branch develop --workdir /tmp/pdf-check \
       serial_dev 2>&1 | tee out.test_serial_dev

    docker exec -it \
       --env OMP_NUM_THREADS=4 \
       --env OMP_SCHEDULE=dynamic \
       --env OMPI_MCA_btl_vader_single_copy_mechanism=none \
       --env MPIEXEC="mpiexec -n 4 --oversubscribe --allow-run-as-root" \
       ${DOCKER_CONTAINER_NAME} \
       pdf-check.sh --branch develop --workdir /tmp/pdf-check \
       parallel 2>&1 | tee out.test_parallel
}


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
if [ -d build ]; then
    rm -rf build
fi
mkdir build
chmod 777 build

docker rm -f ${DOCKER_CONTAINER_NAME} 2>&1 > /dev/null || true
docker run -d --rm \
    --name ${DOCKER_CONTAINER_NAME} \
    -v "${PWD}:/work/ProteinDF" \
    hiracchi/pdf-builder:${PDF_BUILDER_VER}

#docker exec -it ${CONTAINER_NAME} pdf-checkout.sh --branch ${BRANCH} ProteinDF
docker exec -it ${DOCKER_CONTAINER_NAME} pdf-checkout.sh --branch ${PDF_BRIDGE_BRANCH} ProteinDF_bridge
docker exec -it ${DOCKER_CONTAINER_NAME} pdf-checkout.sh --branch ${PDF_PYTOOLS_BRANCH} ProteinDF_pytools

docker exec -it \
    --env MPIEXEC_FLAGS='--oversubscribe;--allow-run-as-root' \
    ${DOCKER_CONTAINER_NAME} pdf-build.sh --srcdir /work/ProteinDF 2>&1 | tee pdf-build.ProteinDF.log
docker exec -it ${DOCKER_CONTAINER_NAME} pdf-build.sh --srcdir /work/ProteinDF_bridge
docker exec -it ${DOCKER_CONTAINER_NAME} pdf-build.sh --srcdir /work/ProteinDF_pytools


run_test

#docker exec -it \
#    ${DOCKER_CONTAINER_NAME} /bin/bash
