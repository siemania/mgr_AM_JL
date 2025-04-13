#!/usr/bin/env bash

mgl_dir="MGLTools"
mgl_tar="${mgl_dir}.tar.gz"

# download MGLTools
wget -O "$mgl_tar" "https://ccsb.scripps.edu/mgltools/download/491/"

# unpack
tar -xvzf "$mgl_tar"
rm "$mgl_tar"

# move directories
mgl_versioned_name="mgltools_x86_64Linux2_1.5.7"
mv "$mgl_versioned_name" "$mgl_dir"

# install MGLTools
cd "$mgl_dir" || exit
source install.sh
alias pythonsh="$mgl_dir/bin/pythonsh"
cd ..

#download wheel
wget https://files.pythonhosted.org/packages/59/b0/11710a598e1e148fb7cbf9220fd2a0b82c98e94efbdecb299cb25e7f0b39/wheel-0.33.6.tar.gz
pythonsh -m easy_install  wheel-0.33.6.tar.gz
rm wheel-0.33.6.tar.gz
#download pip 20.3
wget https://files.pythonhosted.org/packages/03/41/6da553f689d530bc2c337d2c496a40dc9c0fdc6481e5df1f3ee3b8574479/pip-20.3.tar.gz
tar -xzvf pip-20.3.tar.gz
cd pip-20.3/ || exit
pythonsh setup.py install
cd ..
rm -rd pip-20.3

source install_packages.sh