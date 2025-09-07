#!/bin/bash

# rsync -av --dry-run --delete /run/media/enantiomer-s/Analysis/ComparativeGenomics4AL keen:/nfs_share/wyb/ShareFromHuang
# rsync -av --dry-run --delete /run/media/enantiomer-s/Analysis/ComparativeGenomics4AL wyb:/Users/wangyibang/Documents/01research/bianwenyu/ShareFromHuang
rsync -av --dry-run --delete /run/media/enantiomer-s/BacLib/ComparativeGenomics4AL_Publication keen:/nfs_share/huang/SharedProject/

echo "Do you wish to push those files to lab server?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) echo "--- rsync start ---"; rsync -av --delete /run/media/enantiomer-s/BacLib/ComparativeGenomics4AL_Publication keen:/nfs_share/huang/SharedProject/; break;;
        No ) exit;;
    esac
done
