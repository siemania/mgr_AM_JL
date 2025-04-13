#!/usr/bin/env bash

packages=(
  "tqdm 4.36.1"
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