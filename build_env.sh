#!/bin/zsh

# Path to conda installation
CONDA_PATH=$(conda info --base)

# Make conda tools visible
source ${CONDA_PATH}/etc/profile.d/conda.sh

# Environment names 
MINC_ENV_NAME="clustering-autism-minc-env"
ENV_NAME="clustering-autism-env"

# Create MINC environment and install minc-toolkit-v2
echo "Building MINC environment..."
conda create -n $MINC_ENV_NAME -c minc-forge minc-toolkit-v2 -y

# Create conda environment and install R
echo "Building R environment..."
conda create -n $ENV_NAME -c conda-forge r-base=4.4.2 python=3.12.10 -y

# Environment paths
MINC_ENV_PATH=${CONDA_PATH}/envs/${MINC_ENV_NAME}
ENV_PATH=${CONDA_PATH}/envs/${ENV_NAME}

# Set $MINC_TOOLKIT upon activation 
cat <<EOF > ${ENV_PATH}/etc/conda/activate.d/activate-minc-toolkit.sh
if [ -n "\${MINC_TOOLKIT}" ]; then
  export MINC_TOOLKIT_PREV="\${MINC_TOOLKIT}"
fi
export MINC_TOOLKIT=$MINC_ENV_PATH
EOF

# Unset $MINC_TOOLKIT upon deactivation
cat <<EOF > ${ENV_PATH}/etc/conda/deactivate.d/deactivate-minc-toolkit.sh
# restore pre-existing MINC_TOOLKIT 
if [ -n "\${MINC_TOOLKIT_PREV}" ]; then
  export MINC_TOOLKIT="\${MINC_TOOLKIT_PREV}"
  unset MINC_TOOLKIT_PREV
else
  unset MINC_TOOLKIT
fi
EOF


# Activate conda environment
conda activate $ENV_NAME

# Install compiled R packages via conda (faster)
#echo "Installing R packages from conda..."
conda install -c conda-forge --file conda_package_list.txt -y
#conda install -c conda-forge r-devtools=2.4.5 r-tidyverse=2.0.0 r-lme4=1.1_37 r-visnetwork=2.1.2 r-rjson=0.2.23 r-dt=0.33 r-doparallel=1.0.17 r-parallelly=1.43.0 r-cowplot=1.1.3 r-rspectra=0.16.2 r-rcpp=1.0.14 r-rcppeigen=0.3.4.0.2 r-rcpparmadillo=14.4.1_1 -y

# Install Rcpp and others from source
# NOTE: Conda installation seems fine.
#echo "Installing R packages from source..."
#Rscript -e 'install.packages(c("Rcpp", "RcppEigen", "RcppArmadillo"), repos="https://cloud.r-project.org")'

# Install RMINC
#echo "Installing RMINC..."
#Rscript -e 'devtools::install_github("Mouse-Imaging-Centre/RMINC", ref = "57ef9122311d255f24c44571f9c68972c1c3cc4f", upgrade = "never")'


#export PATH="\${MINC_TOOLKIT}/pipeline:\${MINC_TOOLKIT}/bin:\$PATH"
