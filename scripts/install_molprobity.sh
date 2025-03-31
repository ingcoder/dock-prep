#!/bin/bash

# Auto-executable functionality
# Check if script has executable permissions, if not, set them
if [ ! -x "$0" ]; then
    echo "Setting executable permissions for this script..."
    chmod +x "$0"
    echo "Please run the script again."
    exit 0
fi

# Define colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# install_molprobity.sh
# Script to install MolProbity and set up environment variables

echo "======================================================"
echo "Installing MolProbity for protein structure validation"
echo "======================================================"


# Set the installation directory to be within pdbqt-converter
echo -e "${GREEN}=== Installing MolProbity ===${NC}"

echo -e "${YELLOW}Checking for MolProbity installation...${NC}"
INSTALL_DIR="$PWD/MolProbity"
if [ -d "$INSTALL_DIR" ]; then
    echo "MolProbity directory already exists at $INSTALL_DIR"
    echo "Would you like to remove it and reinstall? (y/n)"
    read answer
    if [ "$answer" = "y" ]; then
        echo "Removing existing MolProbity installation..."
        rm -rf "$INSTALL_DIR"
    else
        echo "Installation aborted. Using existing installation."
        exit 0
    fi
fi

echo -e "${YELLOW}Downloading MolProbity...${NC}"
# Move to the pdbqt-converter directory and clone MolProbity repository
echo "Cloning MolProbity repository into $PWD..."
git clone https://github.com/rlabduke/MolProbity.git
cd MolProbity


echo -e "${YELLOW}Installing MolProbity...${NC}"
# Make the installation script executable
echo "Setting up installation script..."
chmod +x ./install_via_bootstrap.sh
# Run the setup script (using method 4 for local installation)
echo "Running MolProbity installation script..."
./install_via_bootstrap.sh 4

# Add MolProbity bin to PATH
# echo "export PATH=\$PATH:$INSTALL_DIR/molprobity/bin" >> ~/.zshrc

# 3. Create environment variables for MGLTools
echo -e "${YELLOW}Setting up environment variables...${NC}"
echo "# MolProbity Environment Setup" > molprobity_env.sh
echo "export MOLPROBITY_ROOT=\"$INSTALL_DIR\"" >> molprobity_env.sh
echo "export PATH=\"$INSTALL_DIR/bin:\$PATH\"" >> molprobity_env.sh

# Make the environment file executable
chmod +x mgltools_env.sh

echo ""
echo "======================================================"
echo "MolProbity installation complete!"
echo "======================================================"
echo ""
echo "MolProbity has been installed in: $INSTALL_DIR"
echo ""
echo "To activate the changes, please run:"
echo "  source $SHELL_CONFIG"
echo ""
echo "You can now use MolProbity tools like 'reduce' and 'probe'"
echo "from the command line."
echo "======================================================" 