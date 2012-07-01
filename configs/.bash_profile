#/Applications/Gnuplot.app/Contents/Resources/bin

ADD=/usr/local/sbin
PATH=$ADD:/Library/Frameworks/GTK+.framework/Resources/bin:$PATH
export PATH=$PATH:"/usr/local/bin:usr/local/libexec/gcc/i386-apple-darwin9.0.0/4.2.3:/usr/local/gromacs/bin:/usr/local/gromacs-3/bin:/Users/grace/Documents/projectIG/android-sdk-mac_x86-1.1_r1/tools:/sw/bin:/Users/grace/projects/cakephp/cake/console"


export GDFONTPATH=$HOME/Library/Fonts:/Library/Fonts:/System/Library/Fonts
export LSCOLORS=ExFxCxDxBxegedabagacad
export DSSP=/Volumes/IomegaHDD/home/grace/src/dssp_src/dssp/dsspcmbi

test -r /sw/bin/init.sh && . /sw/bin/init.sh


alias new="ls -alrt"

# Your previous /Users/grace/.bash_profile file was backed up as /Users/grace/.bash_profile.macports-saved_2010-01-16_at_16:43:25
##
# MacPorts Installer addition on 2010-01-16_at_16:43:25: adding an appropriate PATH variable for use with MacPorts.
#export PATH=/opt/local/bin:/opt/local/sbin:$PATH
# Finished adapting your PATH environment variable for use with MacPorts.
# Setting PATH for Python 2.7
# The orginal version is saved in .bash_profile.pysave
#PATH="/Library/Frameworks/Python.framework/Versions/2.7/bin:${PATH}"
#export PATH


# Setting PATH for MacPython 2.6
# The orginal version is saved in .bash_profile.pysave
PATH="/Library/Frameworks/Python.framework/Versions/2.6/bin:$HOME/bin:${PATH}"
export PATH

# git completion
source ~/bin/git-completion.bash
export TERM="xterm-color"
alias is="is -G"
PS1="\[\e[01;32m\]\u@\h\[\e[01;34m\]\w\$\[\e[0;33;49m\]\$(__git_ps1)\[\033[00m\] "
export CLICOLOR=1

source ~/.bashrc
