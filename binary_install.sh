#!/bin/bash
#version 0.1

#Information in case we can't find the correct binary
help()
{
  echo "Please make sure you are running the latest version of this install script"
  echo "You can always build PHAST from source, it only takes about 5 minutes to do so, just visit http://compgen.bscb.cornell.edu/phast/ for instructions"
  updatetype=none
}


#Make sure the user is root, otherwise we can not change necessary files
if [ `id -u` != '0' ]; then
  echo  "ERROR: Phast requires root privilages to install"
  echo  "ERROR: Please re-run this installer with root privilages"
  exit 1
fi

if [ -z "$1" ]; then
  echo "Making sure you have the latest installer"
  rm -f phast_latest_install.sh
  wget http://compgen.bscb.cornell.edu/phast/install.sh -O phast_latest_install.sh
  if [ $? -eq 0 ]; then
    if [ -f phast_latest_install.sh ]; then
      sh phast_latest_install.sh update
      rm -f phast_latest_install.sh 
      exit $?
    fi
  else
    curl http://compgen.bscb.cornell.edu/phast/install.sh -o phast_latest_install.sh
    if [ $? -eq 0 ]; then
      if [ -f phast_latest_install.sh ]; then
        sh phast_latest_install.sh update
        rm -f phast_latest_install.sh
        exit $?
      fi
    else
      echo "Could not check for installer updates"
    fi
  fi
fi

#If we have previously run this script, remove old phast repos
if [ -f /etc/apt/sources.list ]; then
  sed 's/.*phast.*//g' -i /etc/apt/sources.list
elif [ -f /etc/yum.conf ]; then
  awk '/\[all\]/{c=3}!(c&&c--)' /etc/yum.conf > /etc/yumtemp.conf && mv /etc/yumtemp.conf /etc/yum.conf
fi

#dectect OS and add correct link to compgen repository
if [ -f /etc/apt/sources.list ]; then
  updatetype=apt
  echo "Apt-Get compatable system"
  echo "deb http://compgen.bscb.cornell.edu/phast/apt all free" >> /etc/apt/sources.list
elif [ -f /etc/yum.conf ]; then
  updatetype=yum
  echo "Detected yum compatable system"
  echo "[phast]" >> /etc/yum.conf
  echo "baseurl=http://compgen.bscb.cornell.edu/phast/yum" >> /etc/yum.conf
  echo "gpgcheck=0" >> /etc/yum.conf
elif [ "`cat /etc/*release | grep -i opensuse`" != "" ]; then
  updatetype=zypper
  echo "Detected zypper compatable system"
  zypper sa -t YUM http://compgen.bscb.cornell.edu/phast/yum all
elif [ "`cat /etc/*release | grep -i mandriva`" != "" ]; then
  echo "Running on Mandriva, urpmi repository not yet available"
  help
else
  echo "Sorry, there is no pre-compiled binaries for your distribution."
  help
fi

#Perform the actual install of phast
if [ "$updatetype" = "apt" ]; then
  echo "We need to run 'apt-get update' now"
  echo "Press enter to proceed"
  read proceed
  apt-get update
  echo "Now installing phast via 'apt-get install phast'"
  apt-get install phast
elif [ "$updatetype" = "yum" ]; then
  echo "Now installing phast via 'yum install phast'"
  yum install phast
elif [ "$updatetype" = "zypper" ]; then
  echo "Now installing phast via 'zypper in phast'"
  zypper in phast
fi

