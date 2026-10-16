#!/bin/bash

set -eo pipefail
trap 'echo "Setup failed at line $LINENO" >&2' ERR

# Path to conda installation
CONDA_PATH=$(conda info --base)

# Make conda tools visible
source "${CONDA_PATH}"/etc/profile.d/conda.sh

# Environment names 
MINC_ENV_NAME="minc-env"
ENV_NAME="clustering-autism-env"

# Environment paths
MINC_ENV_PATH=${CONDA_PATH}/envs/${MINC_ENV_NAME}
ENV_PATH=${CONDA_PATH}/envs/${ENV_NAME}

# Create MINC environment and install minc-toolkit-v2
if [ -d "${MINC_ENV_PATH}/conda-meta" ]; then
  echo -e "\nMINC environment ${MINC_ENV_NAME} already exists; skipping creation."
else
  echo -e "\nBuilding MINC environment..."
  conda create -n "$MINC_ENV_NAME" -c minc-forge minc-toolkit-v2 -y
fi

# Validate both newly created and reused MINC environments.
for tool in autocrop perl mincinfo mincresample mincreshape; do
  if [ ! -x "${MINC_ENV_PATH}/bin/${tool}" ]; then
    echo "Missing MINC dependency: ${MINC_ENV_PATH}/bin/${tool}" >&2
    echo "Repair the ${MINC_ENV_NAME} environment before rerunning setup." >&2
    exit 1
  fi
done

(
  export PATH="${MINC_ENV_PATH}/bin:$PATH"
  export PERL5LIB="${MINC_ENV_PATH}/perl:${MINC_ENV_PATH}/pipeline${PERL5LIB:+:$PERL5LIB}"

  autocrop -version
  mincinfo -version
  mincresample -version
  mincreshape -version
)

# Create conda environment and install R
echo -e "\nBuilding R environment..."
conda create -n $ENV_NAME -c conda-forge r-base=4.4.2 python=3.12.10 -y

# Install compiled R packages via conda (faster)
echo -e "\nInstalling R packages..."
conda install -n $ENV_NAME -c conda-forge --file R_packages_conda.txt -y

echo -e "\nInstalling python packages..."
conda install -n $ENV_NAME -c conda-forge --file python_packages_conda.txt -y

# Set $MINC_TOOLKIT upon activation 
cat <<EOF > "${ENV_PATH}"/etc/conda/activate.d/activate-minc-toolkit.sh
if [ -n "\${MINC_TOOLKIT}" ]; then
  export MINC_TOOLKIT_PREV="\${MINC_TOOLKIT}"
fi
export MINC_TOOLKIT=$MINC_ENV_PATH
export CLUSTERING_AUTISM_PERL5LIB_WAS_SET=\${PERL5LIB+x}
export CLUSTERING_AUTISM_PERL5LIB_PREV="\${PERL5LIB-}"
export PERL5LIB="\${MINC_TOOLKIT}/perl:\${MINC_TOOLKIT}/pipeline\${PERL5LIB:+:\$PERL5LIB}"
EOF

# Unset $MINC_TOOLKIT upon deactivation
cat <<EOF > "${ENV_PATH}"/etc/conda/deactivate.d/deactivate-minc-toolkit.sh
# restore pre-existing MINC_TOOLKIT 
if [ -n "\${MINC_TOOLKIT_PREV}" ]; then
  export MINC_TOOLKIT="\${MINC_TOOLKIT_PREV}"
  unset MINC_TOOLKIT_PREV
else
  unset MINC_TOOLKIT
fi
if [ "\${CLUSTERING_AUTISM_PERL5LIB_WAS_SET}" = x ]; then
  export PERL5LIB="\${CLUSTERING_AUTISM_PERL5LIB_PREV}"
else
  unset PERL5LIB
fi
unset CLUSTERING_AUTISM_PERL5LIB_WAS_SET
unset CLUSTERING_AUTISM_PERL5LIB_PREV
EOF

# Activate conda environment
conda activate $ENV_NAME

# Install RMINC from Github
echo -e "\nInstalling RMINC from Github..."
Rscript -e 'devtools::install_github("Mouse-Imaging-Centre/RMINC", ref = "57ef9122311d255f24c44571f9c68972c1c3cc4f", upgrade = "never")'

# Install MRIcrotome from Github
echo -e "\nInstalling MRIcrotome from Github..."
Rscript -e 'devtools::install_github("Mouse-Imaging-Centre/MRIcrotome", upgrade = "never")'

# Install Python datatable and pyminc using pip 
echo -e "\nInstalling python packages using pip..."
pip3 install -r python_packages_pip.txt

# Deactivate conda environment
conda deactivate

# Build hook to set environment paths upon activation
cat <<EOF > "${ENV_PATH}"/etc/conda/activate.d/set-env-vars.sh
if [ -n "\${PYTHONPATH}" ]; then
  export PYTHONPATH_PREV="\${PYTHONPATH}"
fi
export PROJECTPATH="$PWD"
export SRCPATH="\$PROJECTPATH/src"
export PYTHONPATH="\${SRCPATH}"
export PATH_PREFIX="\${SRCPATH}:\${SRCPATH}/drivers:\${SRCPATH}/pipelines:\${MINC_TOOLKIT}/pipeline:\${MINC_TOOLKIT}/bin"
export PATH="\${CONDA_PREFIX}/bin:\${PATH_PREFIX}:\$PATH"
EOF

# Build hook to unset environment paths upon deactivation
cat <<EOF > "${ENV_PATH}"/etc/conda/deactivate.d/unset-env-vars.sh
if [ -n "\${PYTHONPATH_PREV}" ]; then
  export PYTHONPATH="\${PYTHONPATH_PREV}"
  unset PYTHONPATH_PREV
else
  unset PYTHONPATH
fi
export PATH=\${PATH#"\${PATH_PREFIX}:"}
unset PROJECTPATH
unset SRCPATH
unset PATH_PREFIX
EOF
