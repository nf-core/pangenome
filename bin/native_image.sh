#!/bin/bash

if [ ! -d pggb ]; then
    git clone git@github.com:pangenome/pggb.git
fi

if [ ! -d pangenome ]; then
    git clone git@github.com:nf-core/pangenome.git
fi

if [ ! -d pangenome_image ]; then
    mkdir pangenome_image
fi

cd pangenome_image
cp ../pggb/Dockerfile .
# we compile everything in Release mode
sed -i "s/-DCMAKE_BUILD_TYPE=Generic/-DCMAKE_BUILD_TYPE=Release/g" Dockerfile
sed -i "s/-DCMAKE_BUILD_TYPE=Debug/-DCMAKE_BUILD_TYPE=Release/g" Dockerfile
sed -i "s/-DEXTRA_FLAGS='-march=sandybridge -O2'//g" Dockerfile
sed -i "49 s/./#&/" Dockerfile
sed -i "70 s/./#&/" Dockerfile
sed -i "71 s/./#&/" Dockerfile
sed -i "72 s/./#&/" Dockerfile
cp ../pggb/pggb .
if [ ! -d .git ]; then
    cp -r ../pggb/.git .
fi

echo "Kuhl!"

#### SKIP Comment the following lines out if you only want to build an image for PGGB. ####
cat ../pangenome/Dockerfile | tail -n +5 >> Dockerfile
if [ ! -d bin ]; then
    mkdir bin
fi
cp ../pangenome/bin/split_approx_mappings_in_chunks.py bin/
cp ../pangenome/bin/paf2net.py bin/
cp ../pangenome/bin/net2communities.py bin/
cp ../pangenome/environment.yml .
#### SKIP ####

docker build -t ${USER}/pangenome-dev:latest .
#### SKIP Comment the following lines out if you only want to build an image for docker. ####
docker run -d -p 5000:5000 --name registry registry:2
docker image tag $(docker images | grep pangenome-dev | grep latest | grep ${USER} | tr -s ' ' | cut -f 3 -d ' ') localhost:5000/pangenome-dev
docker push localhost:5000/pangenome-dev
SINGULARITY_NOHTTPS=true singularity build pangenome-dev.img docker://localhost:5000/pangenome-dev
docker container stop registry && docker container rm -v registry
#### SKIP ####