#!/bin/bash
#version 0.1

#Information in case we can't find the correct binary
help()
{
  echo "Please make sure you are running the latest version of this install script"
  echo "If your distro has a package called 'lsb' (Linux Standard Base) version 4 or later, you may install it and use our distro independent binaries available at http://compgen.bscb.cornell.edu/phast/repositories/LSB4/"
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
  awk '/\[phast\]/{c=3}!(c&&c--)' /etc/yum.conf > /etc/yumtemp.conf && mv /etc/yumtemp.conf /etc/yum.conf
fi

#dectect OS and add correct link to compgen repository
if [ -f /etc/debian_version ] && [ ! -f /etc/lsb-release ]; then
  updatetype=apt
  echo "Detected Debian"
  if [ "`cat /etc/debian_version | grep ^5\.`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/Debian_5.0  /" >> /etc/apt/sources.list
  elif [ "`cat /etc/debian_version | grep ^4\.`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/Debian_Etch  /" >> /etc/apt/sources.list
  else
    echo "There is no pre-compiled binary for your version of Debian"
    help
  fi
elif [ "`cat /etc/*release | grep -i buntu`" != "" ]; then
  updatetype=apt
  echo "Detected Ubuntu"
  if [ "`cat /etc/*release | grep 10.04`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/xUbuntu_10.04  /" >> /etc/apt/sources.list
  elif [ "`cat /etc/*release | grep 9.10`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/xUbuntu_9.10  /" >> /etc/apt/sources.list
  elif [ "`cat /etc/*release | grep 9.04`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/xUbuntu_9.04  /" >> /etc/apt/sources.list
  elif [ "`cat /etc/*release | grep 8.04`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/xUbuntu_8.04  /" >> /etc/apt/sources.list
  else
    echo "There is no pre-compiled binary for your version of Ubuntu"
    help
  fi
elif [ "`cat /etc/*release | grep -i mint`" != "" ]; then
  updatetype=apt
  echo "Detected Mint"
  if [ "`cat /etc/*release | grep 9`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/xUbuntu_10.04  /" >> /etc/apt/sources.list
  elif [ "`cat /etc/*release | grep 8`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/xUbuntu_9.10  /" >> /etc/apt/sources.list
  elif [ "`cat /etc/*release | grep 7`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/xUbuntu_9.04  /" >> /etc/apt/sources.list
  elif [ "`cat /etc/*release | grep 5`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/xUbuntu_8.04  /" >> /etc/apt/sources.list
  elif [ "`cat /etc/*release | grep 'DISTRIB_RELEASE=1'`" != "" ]; then
    echo "deb http://compgen.bscb.cornell.edu/phast/repositories/Debian_5.0  /" >> /etc/apt/sources.list
  else
    echo "There is no pre-compiled binary for your version of Mint"
    help
  fi

elif [ "`cat /etc/*release | grep -i fedora`" != "" ]; then
  updatetype=yum
  echo "Detected Fedora"
  if [ "`cat /etc/*release | grep 13`" != "" ]; then
    echo "[phast]" >> /etc/yum.conf
    echo "baseurl=http://compgen.bscb.cornell.edu/phast/repositories/Fedora_13" >> /etc/yum.conf
    echo "gpgcheck=0" >> /etc/yum.conf
  elif [ "`cat /etc/*release | grep 12`" != "" ]; then
    echo "[phast]" >> /etc/yum.conf
    echo "baseurl=http://compgen.bscb.cornell.edu/phast/repositories/Fedora_12" >> /etc/yum.conf
    echo "gpgcheck=0" >> /etc/yum.conf
  else
    echo "There is no pre-compiled binary for your version of Fedora"
    help
  fi
elif [ "`cat /etc/*release | grep -i opensuse`" != "" ]; then
  updatetype=zypper
  echo "Detected openSUSE"
  if [ "`cat /etc/*release | grep 11.3`" != "" ]; then
    zypper sa -t YUM http://compgen.bscb.cornell.edu/phast/repositories/openSUSE_11.3/ /
  elif [ "`cat /etc/*release | grep 11.2`" != "" ]; then
    zypper sa -t YUM http://compgen.bscb.cornell.edu/phast/repositories/openSUSE_11.2/ /
  elif [ "`cat /etc/*release | grep 11.1`" != "" ]; then
    zypper sa -t YUM http://compgen.bscb.cornell.edu/phast/repositories/openSUSE_11.1/ /
  else
    echo "There is no pre-compiled binary for your version of openSUSE"
    help
  fi
elif [ "`cat /etc/*release | grep -i mandriva`" != "" ]; then
  echo "Running on Mandriva, urpmi repository not yet available"
elif [ "`cat /etc/*release | grep -i centos`" != "" ]; then
  updatetype=yum
  echo "Detected CentOS"
  if [ "`cat /etc/*release | grep 5\.`" != "" ]; then
    echo "[phast]" >> /etc/yum.conf
    echo "baseurl=http://compgen.bscb.cornell.edu/phast/repositories/CentOS_CentOS-5" >> /etc/yum.conf
    echo "gpgcheck=0" >> /etc/yum.conf
  else
    echo "There is no pre-compiled binary for your version of CentOS"
    help
  fi
elif [ "`cat /etc/*release | grep -i 'red hat'`" != "" ]; then
  updatetype=yum
  echo "Running on Red Hat"
  if [ "`cat /etc/*release | grep 5\.`" != "" ]; then
    echo "[phast]" >> /etc/yum.conf
    echo "baseurl=http://compgen.bscb.cornell.edu/phast/repositories/RedHat_RHEL-5" >> /etc/yum.conf
    echo "gpgcheck=0" >> /etc/yum.conf
  elif [ "`cat /etc/*release | grep 4\.`" != "" ]; then
    if [ -f /etc/yum.conf ]; then
      echo "[phast]" >> /etc/yum.conf
      echo "baseurl=http://compgen.bscb.cornell.edu/phast/repositories/RedHat_RHEL-4" >> /etc/yum.conf
      echo "gpgcheck=0" >> /etc/yum.conf
    else
      echo "RHEL4 does not have yum installed, either install yum and re-run this script, or install via RPM"
      echo "RPM Link: http://compgen.bscb.cornell.edu/phast/repositories/RedHat_RHEL-4"
      help
    fi
  else
    echo "There is no pre-compiled binary for your version of Red Hat"
    help
  fi
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

