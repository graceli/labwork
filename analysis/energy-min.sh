# $1 -- coordinate file
# $2 -- topology file

grompp -f /home/grace/mdp/em.mdp -c $1 -p $2 -o em.tpr
mdrun -s em.tpr -deffnm em

