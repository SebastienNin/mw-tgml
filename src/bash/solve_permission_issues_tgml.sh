# Script to execute by all users who have files in this directory.
echo "Fixing permissions for $USER"
# We want rx for other for subdirectories because else Salva sshfs mount point on his Mac can not access these directories.
find /gpfs/tgml -user $USER -type d -exec chmod 775 {} \; -exec chgrp tgml {} \; -exec chmod g+s {} \; -exec setfacl -d -m g::rwx {} \; -exec setfacl -d -m u::rwx {} \; -exec setfacl -d -m o::rx {} \; 
find /gpfs/tgml -user $USER -type f -exec chmod u+rw,g+rw,o+r {} \; -exec chgrp tgml {} \;
