#!/bin/bash

export PYTHONPATH=/Users/grace/projects/exman:/Users/grace/labwork/lib
PYTHON27=/Library/Frameworks/Python.framework/Versions/2.7/bin/python
export PATH=/Users/grace/projects/exman/exman/bin:/Users/grace/Library/Haskell/ghc-7.0.3/lib/pandoc-1.8.1.1/bin:$PYTHON27:$PATH

export GMXLIB=/usr/local/gromacs/share/gromacs/top

SSH_ENV="$HOME/.ssh/environment"

function load_rails_rvm {
    [[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"  # This loads RVM into a shell session.
    rvm use 1.9.2@rails3 > /dev/null
}

function start_agent {
  echo "Initializing new SSH agent..."
  ssh-agent | sed 's/^echo/#echo/' > "${SSH_ENV}"
  echo succeeded
  chmod 600 "${SSH_ENV}"
  . "${SSH_ENV}" > /dev/null
  ssh-add ~/.ssh/github_graceli
  ssh-add ~/.ssh/github_glacier
  ssh-add ~/.ssh/cluster
  ssh-add ~/.ssh/cluster2

}

# Source my aliases
if [ -e "$HOME/.bash_aliases" ]; then
    . ~/.bash_aliases
fi


# Source SSH settings, if applicable
if [ -f "${SSH_ENV}" ]; then
  . "${SSH_ENV}" > /dev/null
  #ps ${SSH_AGENT_PID} doesn't work under cywgin
  ps -ef | grep ${SSH_AGENT_PID} | grep ssh-agent$ > /dev/null || {
    start_agent;
  }
else
  start_agent;
fi

load_rails_rvm

