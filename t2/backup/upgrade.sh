tar cvfz etc_backup.tgz /etc
tar cvfz var_dpkg.tgz /var/lib/dpkg
tar cvfz pkgstates.tgz /var/lib/aptitude/pkgstates
dpkg --get-selections "*" > dpkg_get_selections.txt
