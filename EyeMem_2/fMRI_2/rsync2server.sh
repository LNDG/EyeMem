#!/bin/bash

#tardis_dir="$HOME/eyemem_tardis/data/mri/"
#server_dir="$HOME/eyemem/data/mri"

tardis_dir="$HOME/eyemem_tardis/study_information/A_scripts/G_Git/"
server_dir="$HOME/eyemem/study_information/A_scripts/G_Git"

rsync -rv $tardis_dir $server_dir