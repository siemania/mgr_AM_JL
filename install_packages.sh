#!/usr/bin/env bash

packages=(
  "tqdm 4.36.1"
  "pathlib 1.0.1"
)


# Function to download and install a package
download_and_install() {
  local name="$1"
  local version="$2"
  echo "Processing $name==$version"

  tar_name="${name}-${version}.tar.gz"
  download_url="https://files.pythonhosted.org/packages/source/${name:0:1}/${name}/${tar_name}"

  wget "$download_url"

  pythonsh -m pip install "$tar_name"

  rm "$tar_name"
}

for entry in "${packages[@]}"; do
  read -r name version <<< "$entry"
  download_and_install "$name" "$version"
done

# wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-ubuntu2404.pinsudo mv cuda-ubuntu2404.pin /etc/apt/preferences.d/cuda-repository-pin-600wget https://developer.download.nvidia.com/compute/cuda/12.9.0/local_installers/cuda-repo-ubuntu2404-12-9-local_12.9.0-575.51.03-1_amd64.debsudo dpkg -i cuda-repo-ubuntu2404-12-9-local_12.9.0-575.51.03-1_amd64.debsudo cp /var/cuda-repo-ubuntu2404-12-9-local/cuda-*-keyring.gpg /usr/share/keyrings/sudo apt-get updatesudo apt-get -y install cuda-toolkit-12-9