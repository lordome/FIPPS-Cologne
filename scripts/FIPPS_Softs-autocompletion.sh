#!/bin/bash

LSTClovBuilder_APP="LSTClovBuilder"

_LSTClovBuilder_completion()
{
  cur=${COMP_WORDS[COMP_CWORD]}

  # List of options
  options=""
  options="$options -conf" 
  options="$options -ow" 
  options="$options -of" 
  options="$options -print"
  options="$options -batch"
  options="$options -nc" 
  options="$options -run" 
  options="$options -help"
  options="$options -nloops"

  local prev
  prev=${COMP_WORDS[COMP_CWORD-1]}

  if [[ "$prev" == "-of" || "$prev" == "-conf" ]]
  then
        COMPREPLY=( $( compgen -W "${files}" -- ${cur}) )
  elif [[ "$prev" == "-run" || "$prev" == "-nloops" ]]
  then
        COMPREPLY=""
  else
        COMPREPLY=( $( compgen -W "$options" -- $cur ) )
  fi
}

complete -o default -o nospace -F _LSTClovBuilder_completion $LSTClovBuilder_APP

EventBuilder_APP="EventBuilder"

_EventBuilder_completion()
{
  cur=${COMP_WORDS[COMP_CWORD]}

  # List of options
  options=""
  options="$options -conf" 
  options="$options -evts" 
  options="$options -ow" 
  options="$options -of" 
  options="$options -print"
  options="$options -batch"
  options="$options -run" 
  options="$options -help" 

  local prev
  prev=${COMP_WORDS[COMP_CWORD-1]}

  if [[ "$prev" == "-of" || "$prev" == "-conf" ]]
  then
        COMPREPLY=( $( compgen -W "${files}" -- ${cur}) )
  elif [[ "$prev" == "-run" || "$prev" == "-evts" ]]
  then
        COMPREPLY=""
  else
        COMPREPLY=( $( compgen -W "$options" -- $cur ) )
  fi
}

complete -o default -o nospace -F _EventBuilder_completion $EventBuilder_APP

G4EventBuilder_APP="G4EventBuilder"

_G4EventBuilder_completion()
{
  cur=${COMP_WORDS[COMP_CWORD]}

  # List of options
  options=""
  options="$options -conf"
  options="$options -evts"
  options="$options -ow"
  options="$options -of"
  options="$options -print"
  options="$options -batch"
  options="$options -run"
  options="$options -help"

  local prev
  prev=${COMP_WORDS[COMP_CWORD-1]}

  if [[ "$prev" == "-of" || "$prev" == "-conf" ]]
  then
        COMPREPLY=( $( compgen -W "${files}" -- ${cur}) )
  elif [[ "$prev" == "-run" || "$prev" == "-evts" ]]
  then
        COMPREPLY=""
  else
        COMPREPLY=( $( compgen -W "$options" -- $cur ) )
  fi
}

complete -o default -o nospace -F _G4EventBuilder_completion $G4EventBuilder_APP

PlotResults_APP="PlotResults"

_PlotResults_completion()
{
  cur=${COMP_WORDS[COMP_CWORD]}

  # List of options
  options=""
  options="$options -conf" 
  options="$options -evts" 
  options="$options -workers" 
  options="$options -ow" 
  options="$options -of" 
  options="$options -print"
  options="$options -run" 
  options="$options -help" 
  options="$options -data-check"

  local prev
  prev=${COMP_WORDS[COMP_CWORD-1]}

  if [[ "$prev" == "-of" || "$prev" == "-conf" ]]
  then
        COMPREPLY=( $( compgen -W "${files}" -- ${cur}) )
  elif [[ "$prev" == "-run" || "$prev" == "-evts" || "$prev" == "-workers" ]]
  then
        COMPREPLY=""
  else
        COMPREPLY=( $( compgen -W "$options" -- $cur ) )
  fi
}

complete -o default -o nospace -F _PlotResults_completion $PlotResults_APP

AnaRaw_APP="AnaRaw"

_AnaRaw_completion()
{
  cur=${COMP_WORDS[COMP_CWORD]}

  # List of options
  options=""
  options="$options -conf"
  options="$options -calibconf"
  options="$options -calib"
  options="$options -file"
  options="$options -id"
  options="$options -init"
  options="$options -track"
  options="$options -qtot"
  options="$options -ow"
  options="$options -help"
  options="$options -run"
  options="$options -time-check"

  local prev
  prev=${COMP_WORDS[COMP_CWORD-1]}

  if [[ "$prev" == "-file" || "$prev" == "-conf" ]]
  then
        COMPREPLY=( $( compgen -W "${files}" -- ${cur}) )
  elif [[ "$prev" == "-run" || "$prev" == "-calibconf" ]]
  then
        COMPREPLY=""
  else
        COMPREPLY=( $( compgen -W "$options" -- $cur ) )
  fi
}

complete -o default -o nospace -F _AnaRaw_completion $AnaRaw_APP

xTalk_APP="xTalk"

_xTalk_completion()
{
  cur=${COMP_WORDS[COMP_CWORD]}

  # List of options
  options=""
  options="$options -conf"
  options="$options -run"
  options="$options -ana"
  options="$options -read"
  options="$options -evts"
  options="$options -of"
  options="$options -ow"
  options="$options -help"
  options="$options -file"
  options="$options -mix"
  options="$options -plot"

  local prev
  prev=${COMP_WORDS[COMP_CWORD-1]}

  if [[ "$prev" == "-file" || "$prev" == "-of" || "$prev" == "-conf" ]]
  then
        COMPREPLY=( $( compgen -W "${files}" -- ${cur}) )
  elif [[ "$prev" == "-evts" || "$prev" == "-run"  ]]
  then
        COMPREPLY=""
  else
        COMPREPLY=( $( compgen -W "$options" -- $cur ) )
  fi
}

complete -o default -o nospace -F _xTalk_completion $xTalk_APP

TkT2ROOT_APP="TkT2ROOT"

_TkT2ROOT_completion()
{
  cur=${COMP_WORDS[COMP_CWORD]}

  # List of options
  options=""
  options="$options -file"
  options="$options -gain"

  local prev
  prev=${COMP_WORDS[COMP_CWORD-1]}

  if [[ "$prev" == "-file" ]]
  then
        COMPREPLY=( $( compgen -W "${files}" -- ${cur}) )
  elif [[ "$prev" == "-gain" ]]
  then
        COMPREPLY=""
  else
        COMPREPLY=( $( compgen -W "$options" -- $cur ) )
  fi
}

complete -o default -o nospace -F _TkT2ROOT_completion $TkT2ROOT_APP
