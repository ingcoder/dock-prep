#!/bin/bash
# Simplified MGLTools Installer Script

# Define colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Installing MGLTools (AutoDock Tools) ===${NC}"

# Detect OS
OS_TYPE="unknown"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS_TYPE="linux"
    echo -e "${YELLOW}Detected Linux operating system${NC}"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS_TYPE="mac"
    echo -e "${YELLOW}Detected macOS operating system${NC}"
else
    echo -e "${RED}Unsupported operating system: $OSTYPE${NC}"
    echo -e "${RED}This script supports Linux and macOS only.${NC}"
    exit 1
fi

# Define installation paths based on OS
CONDA_ENV_NAME="mgltools_py27"

if [[ "$OS_TYPE" == "linux" ]]; then
    DOWNLOAD_URL="https://ccsb.scripps.edu/download/532/"
    DOWNLOAD_FILE="mgltools_x86_64Linux2_1.5.7.tar.gz"
    INSTALL_DIR="$PWD/mgltools_x86_64Linux2_1.5.7"
elif [[ "$OS_TYPE" == "mac" ]]; then
    DOWNLOAD_URL="https://ccsb.scripps.edu/download/529/"
    DOWNLOAD_FILE="mgltools_1.5.7_MacOS-X.tar.gz"
    INSTALL_DIR="$PWD/mgltools_1.5.7_MacOS-X"
fi

# Check for conda installation
echo -e "${YELLOW}Checking for conda installation...${NC}"
if ! command -v conda &> /dev/null; then
    echo -e "${RED}Conda is not installed or not in PATH. Please install Miniconda or Anaconda.${NC}"
    exit 1
fi

# Create or activate conda environment with Python 2.7
echo -e "${YELLOW}Setting up conda environment for MGLTools...${NC}"
if conda env list | grep -q "${CONDA_ENV_NAME}"; then
    echo "Conda environment '${CONDA_ENV_NAME}' already exists. Activating it..."
else
    echo "Creating new conda environment '${CONDA_ENV_NAME}' with Python 2.7..."
    conda create -y -n ${CONDA_ENV_NAME} python=2.7 numpy
    if [ $? -ne 0 ]; then
        echo -e "${RED}Failed to create conda environment.${NC}"
        exit 1
    fi
fi

# Activate the conda environment
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV_NAME}
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to activate conda environment.${NC}"
    exit 1
fi
echo -e "${GREEN}Successfully activated conda environment: ${CONDA_ENV_NAME}${NC}"

# 1. Download MGLTools
echo -e "${YELLOW}Downloading MGLTools...${NC}"
if [ -f "$DOWNLOAD_FILE" ]; then
    echo "MGLTools archive already exists. Removing and re-downloading..."
    rm "$DOWNLOAD_FILE"
fi

echo "Attempting to download from: $DOWNLOAD_URL"
curl -L "$DOWNLOAD_URL" -o "$DOWNLOAD_FILE"
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to download MGLTools.${NC}"
    exit 1
fi 

# 2. Extract the archive
echo -e "${YELLOW}Extracting MGLTools...${NC}"
if [ -d "$INSTALL_DIR" ]; then
    echo "MGLTools directory already exists. Removing it before extraction..."
    rm -rf "$INSTALL_DIR"
fi

tar -xvzf "$DOWNLOAD_FILE"
rm "$DOWNLOAD_FILE"
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to extract MGLTools.${NC}"
    exit 1
fi

# Run the MGLTools installer script
echo -e "${YELLOW}Running MGLTools installer script...${NC}"
cd "$INSTALL_DIR"
chmod +x install.sh
./install.sh
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to install MGLTools.${NC}"
    cd - > /dev/null
    exit 1
fi
cd - > /dev/null

# 3. Create environment variables for MGLTools
echo -e "${YELLOW}Creating/updating MGLTools environment variables...${NC}"

CONFIG_FILE="config_env.json"

# Create or overwrite the config file with MGLTools environment
cat > "$CONFIG_FILE" << EOF
{
    "MGL_ROOT": "$INSTALL_DIR",
    "MGL_PACKAGES": "$INSTALL_DIR/MGLToolsPckgs",
    "MGL_BIN": "$INSTALL_DIR/bin",
    "MGL_ENV_NAME": "$CONDA_ENV_NAME"
}
EOF


# Check if file exists and is not empty
if [ ! -s "$CONFIG_FILE" ]; then
# Create or overwrite the config file with MGLTools environment
    cat > "$CONFIG_FILE" << EOF
    {
        "MGL_ROOT": "$INSTALL_DIR",
        "MGL_PACKAGES": "$INSTALL_DIR/MGLToolsPckgs",
        "MGL_BIN": "$INSTALL_DIR/bin",
        "MGL_ENV_NAME": "$CONDA_ENV_NAME"
    }
EOF
else
    # File exists with content, need to update it
    echo "Processing existing configuration file at $CONFIG_FILE"
    # Using a Python script to properly update JSON with a heredoc for better readability
    python <<EOF
    import json

    # Load existing config
    with open('$CONFIG_FILE', 'r') as f:
        try:
            config = json.load(f)
        except json.JSONDecodeError:
            print('Warning: Invalid JSON in config file. Creating new config.')
            config = {}

    # Update MolProbity settings
    config['MGL_ROOT'] = '$INSTALL_DIR'
    config['MGL_PACKAGES'] = '$INSTALL_DIR/MGLToolsPckgs'
    config['MGL_BIN'] = '$INSTALL_DIR/bin'
    config['MGL_ENV_NAME'] = '$CONDA_ENV_NAME'

    # Write updated config
    with open('$CONFIG_FILE', 'w') as f:
        json.dump(config, f, indent=4)
EOF
    echo "Updated existing configuration file with MolProbity settings"
fi

echo -e "${GREEN}Created/updated ${CONFIG_FILE}${NC}"

echo "# MGLTools Environment Setup" > mgltools_env.sh
echo "# First activate the conda environment" >> mgltools_env.sh
echo "eval \"\$(conda shell.bash hook)\"" >> mgltools_env.sh
echo "conda activate ${CONDA_ENV_NAME}" >> mgltools_env.sh
echo "# Then set MGLTools environment variables" >> mgltools_env.sh
echo "export MGL_ROOT=\"$INSTALL_DIR\"" >> mgltools_env.sh
echo "export PYTHONPATH=\"$INSTALL_DIR/MGLToolsPckgs:\$PYTHONPATH\"" >> mgltools_env.sh
echo "export PATH=\"$INSTALL_DIR/bin:\$PATH\"" >> mgltools_env.sh

# Make the environment file executable
chmod +x mgltools_env.sh

echo -e "${GREEN}MGLTools installation complete in conda environment: ${CONDA_ENV_NAME}!${NC}"
echo -e "To use MGLTools, run: source $PWD/mgltools_env.sh"
echo -e "This will activate the conda environment and set up MGLTools paths."
echo -e "Then you can convert PDB files with:"
echo -e "  $INSTALL_DIR/bin/pythonsh $INSTALL_DIR/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r protein.pdb -o protein.pdbqt"