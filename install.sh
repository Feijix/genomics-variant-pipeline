#!/bin/bash

# ============================================================
# Quick Installer - Copy and paste into your Fedora WSL terminal
# ============================================================
# This downloads the setup script and runs it in one command
# ============================================================

set -e

echo "🧬 Genomics Pipeline Quick Installer for Fedora WSL"
echo "===================================================="
echo ""

# Check we're on Fedora
if [ ! -f /etc/fedora-release ]; then
    echo "❌ This script is designed for Fedora only"
    exit 1
fi

# Check we're not root
if [ "$EUID" -eq 0 ]; then
    echo "❌ Please don't run as root. Use: ./install.sh"
    exit 1
fi

# Install essential tools first
echo "📦 Installing essential tools..."
sudo dnf install -y curl wget git java-17-openjdk

# Install Docker
echo "🐳 Installing Docker..."
sudo dnf install -y dnf-plugins-core
sudo dnf config-manager --add-repo https://download.docker.com/linux/fedora/docker-ce.repo
sudo dnf install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo groupadd docker 2>/dev/null || true
sudo usermod -aG docker $USER

# Install Nextflow
echo "🚀 Installing Nextflow..."
mkdir -p "$HOME/.local/bin"
export PATH="$HOME/.local/bin:$PATH"
curl -s https://get.nextflow.io | bash
mv nextflow "$HOME/.local/bin/"
chmod +x "$HOME/.local/bin/nextflow"

# Add to bashrc
if ! grep -q ".local/bin" "$HOME/.bashrc" 2>/dev/null; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
fi

# Start Docker
echo "🔧 Starting Docker..."
sudo dockerd > /dev/null 2>&1 &
sleep 5

# Pull containers
echo "📦 Pulling bioinformatics containers..."
docker pull quay.io/biocontainers/bwa-mem2:2.2.1--h43eeafb_2 &
docker pull google/deepvariant:1.5.0 &
docker pull quay.io/biocontainers/bcftools:1.20--h8b25389_0 &
docker pull quay.io/biocontainers/samtools:1.20--h50ea8bc_0 &
wait

# Verify
echo ""
echo "✅ Setup Complete!"
echo ""
echo "Verifying installations:"
java -version 2>&1 | head -1
nextflow -version 2>&1 | head -1
docker --version
echo ""

echo "================================================"
echo "🎉 You're ready to run the pipeline!"
echo "================================================"
echo ""
echo "Next steps:"
echo "  1. Restart your terminal (or: source ~/.bashrc)"
echo "  2. Copy the genomics-pipeline folder to ~/ (your home)"
echo "  3. Prepare reference: ./bin/prepare_reference.sh /path/to/chr20.fa"
echo "  4. Update data/samplesheet.csv with your samples"
echo "  5. Run: nextflow run main.nf -profile docker"
echo ""
echo "To start Docker in future sessions:"
echo "  sudo dockerd &"
echo ""